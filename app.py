import csv
from io import StringIO
from itertools import groupby
import re
from statistics import median
from textwrap import fill

import plotly.graph_objects as go
import streamlit as st

from sanger import (
    DNA,
    MAX_READ_BASES,
    MAX_REFERENCE_TEXT,
    assemble_reads,
    base_passes,
    format_alignment,
    map_read,
    mask_read_bases,
    parse_ab1,
    parse_references,
    reference_report,
)


BASE_COLORS = {"A": "#07825C", "C": "#2464BA", "G": "#242D36", "T": "#CB4349"}
READ_LABELS = {"F": "Forward (F)", "R": "Reverse (R)"}
READ_COLORS = {"F": "#0D716C", "R": "#5852A3"}
PLOT_CONFIG = {"displaylogo": False, "modeBarButtonsToRemove": ["lasso2d", "select2d"]}


def draw_quality(reads, threshold):
    rows = []
    figure = go.Figure()
    for direction, read in reads.items():
        quality = read.quality
        if quality is None:
            rows.append({"Chiều": direction, "Base đọc": len(read.sequence),
                         f"Đạt Q{threshold}": "—", "Q trung vị": "—",
                         f"Đoạn Q{threshold} dài nhất": "—"})
            continue
        usable = [base_passes(base, q, threshold) for base, q in zip(read.sequence, quality, strict=True)]
        count = sum(usable)
        longest = max((sum(1 for _ in group) for valid, group in groupby(usable) if valid), default=0)
        rows.append({"Chiều": direction, "Base đọc": len(quality),
                     f"Đạt Q{threshold}": f"{count}/{len(quality)} ({count / len(quality):.1%})",
                     "Q trung vị": f"{median(quality):g}",
                     f"Đoạn Q{threshold} dài nhất": f"{longest} nt"})
        figure.add_scatter(x=list(range(1, len(quality) + 1)), y=list(quality),
                           mode="lines+markers", name=direction,
                           marker={"size": 4}, line={"color": READ_COLORS[direction], "width": 1.4},
                           hovertemplate="Base %{x} · Q%{y}<extra>%{fullData.name}</extra>")
    st.dataframe(rows, hide_index=True)
    if figure.data:
        figure.add_hline(y=threshold, line_dash="dash", line_color="#8B9AA0",
                         annotation_text=f"Q{threshold}", annotation_position="top right")
        figure.update_layout(xaxis_title="Vị trí base trong read gốc", yaxis_title="Điểm Q",
                             height=300, margin={"t": 20, "b": 40, "l": 45, "r": 20},
                             legend={"orientation": "h", "y": 1.08})
        st.plotly_chart(figure, config=PLOT_CONFIG, key="quality-chart")
    if any(read.quality is None for read in reads.values()):
        st.caption("Dấu — chỉ tệp không có điểm Q. Các base này vẫn có thể xuất hiện trong căn chỉnh.")


def show_read_bases(reads, threshold):
    st.caption(f"Mỗi FASTA giữ chiều và vị trí gốc. Base dưới Q{threshold} hoặc chưa xác định được ghi N.")
    filtered_reads = {}
    export_names = {}
    for column, (direction, read) in zip(st.columns(len(reads)), reads.items(), strict=True):
        filtered = mask_read_bases(read, threshold)
        stem = re.sub(r"[^a-zA-Z0-9_-]+", "_", read.name.rsplit(".", 1)[0]).strip("_")[:60] or "read"
        name = f"{stem}_{direction}_Q{threshold}"
        filtered_reads[direction] = filtered
        export_names[direction] = name
        column.download_button(f"Tải FASTA {direction}", f">{name}\n{fill(filtered, width=80)}\n",
                               f"{name}_filtered.fasta", "text/plain", type="primary",
                               on_click="ignore", width="stretch")

    directions = list(reads)
    direction = directions[0]
    if len(directions) > 1:
        direction = st.radio("Xem từng base của", directions, format_func=READ_LABELS.get,
                             horizontal=True, key="base-quality-read")
    read = reads[direction]
    quality = read.quality if read.quality is not None else (None,) * len(read.sequence)
    filtered = filtered_reads[direction]
    unknown = filtered.count("N")
    st.caption(f"{len(filtered)} base · {len(filtered) - unknown} base đạt Q{threshold} · {unknown} vị trí N.")
    if read.quality is None:
        st.info("Tệp thiếu điểm Q. Các vị trí được ghi N trong chuỗi sau lọc.")
    rows = [
        {"Vị trí trong read": position, "Base máy gọi": base, "Điểm Q": q, "Base sau lọc": masked}
        for position, (base, q, masked) in enumerate(zip(read.sequence, quality, filtered, strict=True), start=1)
    ]
    st.dataframe(rows, hide_index=True, height=280)
    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    st.download_button(f"Tải CSV từng base của {direction}", buffer.getvalue().encode("utf-8-sig"),
                       f"{export_names[direction]}_bases.csv", "text/csv", on_click="ignore")


def draw_trace(read, key, focus=None):
    if not read.traces:
        st.info("Tệp này không có đủ bốn kênh chromatogram đã xử lý.")
        return
    figure = go.Figure()
    size = len(read.traces["A"])
    for base, color in BASE_COLORS.items():
        figure.add_scatter(y=list(read.traces[base]), name=base, mode="lines",
                           line={"color": color, "width": 1.3},
                           hovertemplate="Điểm quét %{x} · %{y} RFU<extra>%{fullData.name}</extra>")
    if read.peak_locations is not None:
        points = list(read.peak_locations)
        heights = [max(channel[x] for channel in read.traces.values()) for x in points]
        quality = read.quality if read.quality is not None else ["—"] * len(points)
        figure.add_scatter(
            x=points, y=heights, mode="markers", showlegend=False, name="Base call",
            marker={"color": "#607983", "size": 5, "symbol": "circle-open"},
            customdata=list(zip(range(1, len(points) + 1), read.sequence, quality)),
            hovertemplate="Base %{customdata[0]} · %{customdata[1]} · Q%{customdata[2]}<extra></extra>",
        )
        if focus is not None:
            x = points[focus - 1]
            figure.add_vline(x=x, line_color="#0D716C", line_dash="dot",
                             annotation_text=f"Base {focus}")
            figure.update_xaxes(range=[max(0, x - 75), min(size - 1, x + 75)])
    figure.update_layout(xaxis_title="Điểm quét", yaxis_title="Cường độ (RFU)",
                         margin={"t": 26, "b": 45, "l": 50, "r": 20}, height=320,
                         legend={"orientation": "h", "y": 1.12},
                         uirevision=f"{key}-{focus}-{read.name}")
    st.plotly_chart(figure, key=key, config=PLOT_CONFIG)
    if read.peak_locations is None:
        st.caption("Tệp có tín hiệu RFU nhưng thiếu tọa độ peak hợp lệ để gắn base lên đồ thị.")


def csv_report(positions, mappings, threshold):
    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["type", "ref_position", "after_ref_position", "ref_base",
                     "F_base", "F_original_position", "F_Q",
                     "R_base", "R_original_position", "R_Q", "consensus", "status", "Q_threshold"])
    for pos in positions:
        writer.writerow(["reference", pos.index, "", pos.ref_base,
                         *call_fields(pos.calls.get("F")), *call_fields(pos.calls.get("R")),
                         pos.consensus, pos.status, threshold])
    for direction, mapping in mappings.items():
        last_ref = None
        for column in mapping.columns:
            if column.ref_position is not None:
                last_ref = column.ref_position
                continue
            fields = call_fields(column)
            writer.writerow(["insertion_in_alignment", "", last_ref, "-",
                             *(fields if direction == "F" else ("", "", "")),
                             *(fields if direction == "R" else ("", "", "")),
                             "", f"Chèn trong căn chỉnh {direction}", threshold])
    return buffer.getvalue().encode("utf-8-sig")


def call_fields(call):
    return ("", "", "") if call is None else (call.read_base, call.read_position, call.quality)


def show_alignment(reads, ranked, threshold):
    labels = [f"{i + 1}. {ref.name}" for i, (ref, _) in enumerate(ranked)]
    selected = labels[0]
    if len(labels) > 1:
        selected = st.selectbox("Chuỗi tham chiếu", labels,
                                help="Xếp theo tổng điểm căn chỉnh cục bộ của các read.")
    ref, mappings = ranked[labels.index(selected)]
    if not any(mapping.columns for mapping in mappings.values()):
        st.warning("Không tìm thấy đoạn căn chỉnh đủ dài với chuỗi này.")
        return

    positions = reference_report(ref, mappings, threshold)
    covered = sum(pos.consensus in DNA for pos in positions)
    match = sum(pos.status == "Khớp" for pos in positions)
    different = sum(pos.status == "Khác tham chiếu" for pos in positions)
    left, middle, right = st.columns(3)
    left.metric(f"Có base đạt Q{threshold}", f"{covered}/{len(positions)}")
    middle.metric(f"Khớp ref ở Q{threshold}", f"{match}/{len(positions)}")
    right.metric(f"Khác ref ở Q{threshold}", different)

    for direction, mapping in mappings.items():
        if not mapping.columns:
            st.caption(f"{READ_LABELS[direction]} · Không có đoạn căn chỉnh đủ dài.")
            continue
        orientation = "đảo bổ sung" if mapping.reverse else "chiều gốc"
        identity = mapping.identities / len(mapping.columns)
        st.markdown(f"**{READ_LABELS[direction]}** · Ref {mapping.reference_span} · "
                    f"{mapping.identities}/{len(mapping.columns)} khớp ({identity:.1%}) · {orientation}")
        st.code(format_alignment(mapping), language=None, wrap_lines=False)
    st.caption("| khớp · . khác · - gap · ? base mơ hồ. Tọa độ Read là vị trí trong tệp gốc.")

    with st.expander("Kiểm tra từng vị trí"):
        key = f"ref-position-{labels.index(selected)}"
        index = st.number_input("Vị trí trên tham chiếu", min_value=1,
                                max_value=len(positions), value=1, step=1, key=key)
        row = positions[index - 1]
        calls = []
        for direction in reads:
            base, read_pos, q = call_fields(row.calls.get(direction))
            calls.append({"Chiều": direction, "Base tham chiếu": row.ref_base,
                          "Base quan sát": base or "—",
                          "Điểm Q": str(q) if isinstance(q, int) else "—",
                          "Vị trí trong read": str(read_pos) if isinstance(read_pos, int) else "—"})
        st.dataframe(calls, hide_index=True)
        st.caption(f"Consensus {row.consensus} · {row.status}")
        trace_options = [direction for direction, call in row.calls.items()
                         if call.read_position is not None and reads[direction].peak_locations]
        if trace_options:
            direction = st.radio("Xem peak", trace_options, horizontal=True, key="peak-direction")
            draw_trace(reads[direction], "aligned-peak", focus=row.calls[direction].read_position)

    stem = re.sub(r"[^a-zA-Z0-9_-]+", "_", ref.name).strip("_")[:60] or "reference"
    consensus = "".join(pos.consensus for pos in positions)
    fasta = f">{stem}_ref_coordinates_Q{threshold}\n{fill(consensus, width=80)}\n"
    left, right = st.columns(2)
    left.download_button("Tải bảng CSV", csv_report(positions, mappings, threshold),
                         f"{stem}_comparison.csv", "text/csv", on_click="ignore")
    right.download_button("Tải consensus FASTA", fasta, f"{stem}_consensus.fasta",
                          "text/plain", on_click="ignore")
    st.caption("Consensus theo tọa độ tham chiếu. N là vị trí chưa gọi được hoặc F/R xung đột. "
               "Các base chèn được lưu trong CSV.")


def show_read_consensus(reads, threshold):
    st.caption(f"App lọc hai đầu read theo Q{threshold}, lấy đảo bổ sung của R rồi tìm phần chồng lấn với F.")
    with st.expander("Tiêu chí ghép"):
        left, middle, right = st.columns(3)
        min_overlap = left.number_input("Số cặp base tối thiểu", min_value=1,
                                        max_value=MAX_READ_BASES, value=20)
        min_identity = middle.number_input("Tỷ lệ khớp tối thiểu (%)", min_value=1,
                                            max_value=100, value=85)
        min_matches = right.number_input(f"Số vị trí cả hai đạt Q{threshold} và khớp", min_value=1,
                                         max_value=MAX_READ_BASES, value=8, key="assembly-min-matches")
        st.caption("Tỷ lệ khớp tính trên đoạn chồng lấn, kể cả gap trong đoạn đó.")

    result = assemble_reads(reads["F"], reads["R"], threshold, min_overlap=min_overlap,
                            min_identity=min_identity / 100, min_high_quality_matches=min_matches)
    if result.overlap:
        left, middle, right = st.columns(3)
        left.metric("Cặp base chồng lấn", result.overlap)
        middle.metric("Tỷ lệ khớp", f"{result.identity:.1%}")
        right.metric(f"Cả hai đạt Q{threshold} và khớp", result.high_quality_matches)
        with st.expander("Xem căn chỉnh F/R"):
            blocks = [f"F     {result.forward_aligned[i:i + 54]}\n"
                      f"R RC  {result.reverse_aligned[i:i + 54]}"
                      for i in range(0, len(result.forward_aligned), 54)]
            st.code("\n\n".join(blocks), language=None, wrap_lines=False)
        st.caption("RC là chuỗi đảo bổ sung. Dấu gạch ngang đánh dấu gap hoặc phần ngoài vùng chồng lấn.")

    if result.sequence is None:
        st.warning(result.problem)
        return
    st.subheader("Chuỗi consensus F/R")
    st.caption(f"{len(result.sequence)} nt · {result.sequence.count('N')} vị trí chưa xác định")
    st.code(fill(result.sequence, width=60), language=None)
    fasta = f">FR_consensus_no_reference_Q{threshold}\n{fill(result.sequence, width=80)}\n"
    st.download_button("Tải consensus FASTA", fasta, "FR_consensus_no_reference.fasta",
                       "text/plain", on_click="ignore")


def run_app():
    st.set_page_config(page_title="Phân tích Sanger", layout="centered")
    st.title("Phân tích Sanger")
    st.caption("Tải FASTA riêng cho F và R, kiểm tra chất lượng và xem chromatogram.")

    with st.form("inputs", border=True):
        st.caption("Tải lên một hoặc hai tệp .ab1.")
        left, right = st.columns(2)
        forward = left.file_uploader("Forward (F)", type=["ab1", "abi"])
        reverse = right.file_uploader("Reverse (R)", type=["ab1", "abi"])
        with st.expander("Chuỗi tham chiếu · tùy chọn"):
            reference = st.text_area("Chuỗi DNA hoặc FASTA", height=105,
                                     max_chars=MAX_REFERENCE_TEXT,
                                     placeholder="Dán một chuỗi DNA hoặc nhiều bản ghi FASTA.")
        submitted = st.form_submit_button("Phân tích", type="primary")

    if submitted:
        st.session_state.pop("analysis", None)
        uploads = {d: upload for d, upload in (("F", forward), ("R", reverse)) if upload is not None}
        if not uploads:
            st.info("Hãy chọn ít nhất một tệp .ab1 để phân tích.")
            return

        try:
            reads = {}
            for direction, upload in uploads.items():
                try:
                    reads[direction] = parse_ab1(upload.name, upload.getvalue())
                except ValueError as exc:
                    raise ValueError(f"{READ_LABELS[direction]} · {exc}") from exc
            refs = parse_references(reference)
            ranked = [(ref, {d: map_read(ref, read, d == "R") for d, read in reads.items()}) for ref in refs]
            ranked.sort(key=lambda result: sum(mapping.score for mapping in result[1].values()), reverse=True)
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state.analysis = (reads, ranked)

    if "analysis" not in st.session_state:
        return
    reads, ranked = st.session_state.analysis
    has_pair = "F" in reads and "R" in reads
    threshold_column, _ = st.columns([1, 2])
    threshold = threshold_column.number_input("Ngưỡng chất lượng Q", min_value=1, max_value=60,
                                              value=20, step=1)
    tab_labels = ["Tải FASTA F/R", "Chất lượng", "Chromatogram"]
    if ranked:
        tab_labels.append("So sánh tham chiếu")
    elif has_pair:
        tab_labels.append("Consensus F/R")
    tabs = st.tabs(tab_labels)
    with tabs[0]:
        show_read_bases(reads, threshold)
    with tabs[1]:
        draw_quality(reads, threshold)
    with tabs[2]:
        directions = list(reads)
        direction = st.radio("Chiều đọc", directions, format_func=READ_LABELS.get,
                             horizontal=True) if len(directions) > 1 else directions[0]
        st.text(reads[direction].name)
        draw_trace(reads[direction], "chromatogram")
    if ranked:
        with tabs[3]:
            show_alignment(reads, ranked, threshold)
    elif has_pair:
        with tabs[3]:
            show_read_consensus(reads, threshold)


if __name__ == "__main__":
    run_app()
