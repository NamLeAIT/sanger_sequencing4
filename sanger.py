"""Read Sanger traces and assemble or compare observed DNA bases."""

from dataclasses import dataclass
from functools import lru_cache
from io import BytesIO, StringIO
from struct import error as StructError

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import Array
from Bio.Seq import reverse_complement


DNA = frozenset("ACGT")
IUPAC = DNA | frozenset("RYSWKMBDHVN")
MAX_REFERENCE_TEXT = 60_000
MAX_REFERENCE_BASES = 50_000
MAX_REFERENCE_COUNT = 20
MAX_TRACE_BYTES = 5 * 1024 * 1024
MAX_READ_BASES = 3000


@dataclass(frozen=True)
class Read:
    name: str
    sequence: str
    quality: tuple[int, ...] | None
    peak_locations: tuple[int, ...] | None
    traces: dict[str, tuple[int, ...]]


@dataclass(frozen=True)
class Reference:
    name: str
    sequence: str


@dataclass(frozen=True)
class Column:
    ref_position: int | None
    read_position: int | None
    ref_base: str
    read_base: str
    quality: int | None

    @property
    def symbol(self) -> str:
        if "-" in (self.ref_base, self.read_base):
            return " "
        if self.ref_base not in DNA or self.read_base not in DNA:
            return "?"
        return "|" if self.ref_base == self.read_base else "."


@dataclass(frozen=True)
class Mapping:
    columns: tuple[Column, ...]
    score: float
    reverse: bool

    @property
    def identities(self) -> int:
        return sum(col.symbol == "|" for col in self.columns)

    @property
    def paired_bases(self) -> int:
        return sum(col.ref_position is not None and col.read_position is not None
                   for col in self.columns)

    @property
    def reference_span(self) -> str:
        positions = [col.ref_position for col in self.columns if col.ref_position is not None]
        return f"{positions[0]}–{positions[-1]}" if positions else "—"


@dataclass(frozen=True)
class Position:
    index: int
    ref_base: str
    calls: dict[str, Column]
    consensus: str
    status: str


@dataclass(frozen=True)
class ReadAssembly:
    forward_aligned: str = ""
    reverse_aligned: str = ""
    overlap: int = 0
    identity: float = 0.0
    high_quality_matches: int = 0
    sequence: str | None = None
    problem: str = ""


def parse_ab1(name: str, contents: bytes) -> Read:
    if len(contents) > MAX_TRACE_BYTES:
        raise ValueError("Tệp .ab1 vượt quá 5 MB.")
    if not contents.startswith(b"ABIF") or len(contents) < 34:
        raise ValueError("Tệp không có định dạng ABIF hợp lệ.")
    try:
        record = SeqIO.read(BytesIO(contents), "abi")
    except (ValueError, OSError, KeyError, IndexError, TypeError, StructError) as exc:
        raise ValueError("Không đọc được tệp .ab1. Hãy kiểm tra tệp gốc của máy.") from exc
    sequence = str(record.seq).upper()
    if not sequence or not set(sequence) <= IUPAC:
        raise ValueError("Tệp không chứa chuỗi DNA đã gọi hợp lệ.")
    if len(sequence) > MAX_READ_BASES:
        raise ValueError(f"Read dài hơn giới hạn {MAX_READ_BASES} base.")
    raw = record.annotations.get("abif_raw", {})
    order = raw.get("FWO_1", b"").decode("ascii", errors="ignore")
    channel_data = [raw.get(f"DATA{index}") for index in range(9, 13)]
    traces = {}
    if len(order) == 4 and set(order) == DNA and all(
        isinstance(channel, tuple) and len(channel) >= 2 for channel in channel_data
    ):
        if len({len(channel) for channel in channel_data}) == 1:
            traces = dict(zip(order, channel_data))
    positions = raw.get("PLOC2")
    if isinstance(positions, int):
        positions = (positions,)
    if not (
        traces and isinstance(positions, tuple) and len(positions) == len(sequence)
        and 0 <= positions[0] <= positions[-1] < len(channel_data[0])
        and all(a < b for a, b in zip(positions, positions[1:]))
    ):
        positions = None
    scores = record.letter_annotations.get("phred_quality")
    quality = tuple(scores) if scores is not None and len(scores) == len(sequence) else None
    return Read(name, sequence, quality, positions, traces)


def parse_references(text: str) -> list[Reference]:
    text = text.strip().lstrip("\ufeff")
    if not text:
        return []
    if text.startswith(">"):
        records = [(record.description, str(record.seq)) for record in SeqIO.parse(StringIO(text), "fasta")]
    else:
        records = [("Tham chiếu", text)]
    if not records:
        raise ValueError("Hãy nhập chuỗi DNA hoặc FASTA hợp lệ.")
    if len(records) > MAX_REFERENCE_COUNT:
        raise ValueError(f"Hãy nhập tối đa {MAX_REFERENCE_COUNT} chuỗi tham chiếu.")
    refs = []
    for name, letters in records:
        sequence = "".join(letters.split()).upper()
        if not sequence or set(sequence) - IUPAC:
            raise ValueError("Tham chiếu phải là chuỗi DNA, có thể dùng ký hiệu IUPAC như N.")
        refs.append(Reference(name or f"Tham chiếu {len(refs) + 1}", sequence))
    if sum(len(ref.sequence) for ref in refs) > MAX_REFERENCE_BASES:
        raise ValueError("Tổng độ dài các chuỗi tham chiếu phải từ 50.000 nt trở xuống.")
    return refs


def make_aligner(mode: str) -> PairwiseAligner:
    aligner = PairwiseAligner(mode=mode)
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -2
    # Treat ambiguity codes as neutral during alignment and exclude them from identities.
    aligner.substitution_matrix = ambiguity_matrix()
    return aligner


def map_read(ref: Reference, read: Read, prefer_reverse: bool) -> Mapping:
    aligner = make_aligner("local")
    reverse = reverse_complement(read.sequence)
    forward_score = aligner.score(ref.sequence, read.sequence)
    reverse_score = aligner.score(ref.sequence, reverse)
    is_reverse = reverse_score > forward_score or (reverse_score == forward_score and prefer_reverse)
    query = reverse if is_reverse else read.sequence
    alignments = aligner.align(ref.sequence, query)
    first = next(iter(alignments), None)
    if first is None:
        return Mapping((), 0, is_reverse)
    row_ref, row_read = first.indices
    columns = []
    for ref_index, oriented_index in zip(row_ref, row_read, strict=True):
        r = int(ref_index)
        q = int(oriented_index)
        instrument_index = len(query) - 1 - q if is_reverse else q
        columns.append(Column(
            ref_position=r + 1 if r >= 0 else None,
            read_position=instrument_index + 1 if q >= 0 else None,
            ref_base=ref.sequence[r] if r >= 0 else "-",
            read_base=query[q] if q >= 0 else "-",
            quality=read.quality[instrument_index] if q >= 0 and read.quality is not None else None,
        ))
    mapping = Mapping(tuple(columns), float(first.score), is_reverse)
    # Short incidental alignments are insufficient unless the input itself is short.
    if mapping.paired_bases < min(12, len(ref.sequence), len(query)):
        return Mapping((), 0, is_reverse)
    return mapping


@lru_cache(maxsize=1)
def ambiguity_matrix():
    matrix = Array(alphabet="ACGTRYSWKMBDHVN", dims=2)
    for a in DNA:
        for b in DNA:
            matrix[a, b] = 2 if a == b else -3
    return matrix


def reference_report(ref: Reference, mappings: dict[str, Mapping], threshold: int) -> list[Position]:
    indexed = {
        direction: {col.ref_position: col for col in mapping.columns if col.ref_position is not None}
        for direction, mapping in mappings.items()
    }
    positions = []
    for index, base in enumerate(ref.sequence, start=1):
        calls = {direction: columns[index] for direction, columns in indexed.items() if index in columns}
        high = [col.read_base for col in calls.values()
                if col.read_base in DNA and col.quality is not None and col.quality >= threshold]
        if len(set(high)) > 1:
            consensus, status = "N", "Xung đột F/R"
        elif high:
            consensus = high[0]
            if base not in DNA:
                status = "Ref mơ hồ"
            else:
                status = "Khớp" if consensus == base else "Khác tham chiếu"
        else:
            consensus = "N"
            status = "Q thấp / chưa rõ" if any(col.read_position for col in calls.values()) else "Thiếu base"
        positions.append(Position(index, base, calls, consensus, status))
    return positions


def format_alignment(mapping: Mapping, width: int = 48) -> str:
    blocks = []
    for start in range(0, len(mapping.columns), width):
        columns = mapping.columns[start:start + width]
        ref_pos = [col.ref_position for col in columns if col.ref_position is not None]
        read_pos = [col.read_position for col in columns if col.read_position is not None]
        ref = "".join(col.ref_base for col in columns)
        read = "".join(col.read_base for col in columns)
        middle = "".join(col.symbol for col in columns)
        ref_start, ref_end = (ref_pos[0], ref_pos[-1]) if ref_pos else ("—", "—")
        read_start, read_end = (read_pos[0], read_pos[-1]) if read_pos else ("—", "—")
        blocks.append(f"Ref   {ref_start:>5}  {ref}  {ref_end}\n"
                      f"             {middle}\n"
                      f"Read  {read_start:>5}  {read}  {read_end}")
    return "\n\n".join(blocks)


def base_passes(base: str, quality: int | None, threshold: int) -> bool:
    return base in DNA and quality is not None and quality >= threshold


def mask_read_bases(read: Read, threshold: int) -> str:
    if read.quality is None:
        return "N" * len(read.sequence)
    return "".join(base if base_passes(base, q, threshold) else "N"
                   for base, q in zip(read.sequence, read.quality, strict=True))


def trim_low_quality_ends(read: Read, threshold: int) -> tuple[str, tuple[int, ...]]:
    if read.quality is None:
        return "", ()
    good_indices = [i for i, (base, q) in enumerate(zip(read.sequence, read.quality, strict=True))
                    if base_passes(base, q, threshold)]
    if not good_indices:
        return "", ()
    start, end = good_indices[0], good_indices[-1] + 1
    return read.sequence[start:end], read.quality[start:end]


def assemble_reads(forward: Read, reverse: Read, threshold: int, *,
                   min_overlap: int = 20, min_identity: float = 0.85,
                   min_high_quality_matches: int = 8) -> ReadAssembly:
    if forward.quality is None or reverse.quality is None:
        return ReadAssembly(problem="Cần điểm Q trong cả hai tệp để đánh giá phần chồng lấn.")

    f_sequence, f_quality = trim_low_quality_ends(forward, threshold)
    r_sequence, r_quality = trim_low_quality_ends(reverse, threshold)
    if not f_sequence or not r_sequence:
        return ReadAssembly(problem=f"Ít nhất một read không có base đạt Q{threshold}.")
    r_sequence = reverse_complement(r_sequence)
    r_quality = r_quality[::-1]

    # Permit unpaired ends while keeping penalties for gaps inside the overlap.
    aligner = make_aligner("global")
    aligner.end_gap_score = 0
    alignments = aligner.align(f_sequence, r_sequence)
    alignment = next(alignments)
    f_indices, r_indices = alignment.indices.tolist()
    pairs = [i for i, (f, r) in enumerate(zip(f_indices, r_indices, strict=True)) if f >= 0 and r >= 0]
    if not pairs or alignment.score <= 0:
        return ReadAssembly(problem="Không tìm được phần chồng lấn có điểm căn chỉnh dương.")
    start, end = pairs[0], pairs[-1]
    overlap = len(pairs)
    f_aligned = "".join(f_sequence[i] if i >= 0 else "-" for i in f_indices)
    r_aligned = "".join(r_sequence[i] if i >= 0 else "-" for i in r_indices)
    identity = sum(f_aligned[i] == r_aligned[i] and f_aligned[i] in DNA
                   for i in pairs) / (end - start + 1)

    consensus = []
    high_quality_matches = 0
    for column, (f, r) in enumerate(zip(f_indices, r_indices, strict=True)):
        f_base, r_base = f_aligned[column], r_aligned[column]
        f_high = f >= 0 and f_base in DNA and f_quality[f] >= threshold
        r_high = r >= 0 and r_base in DNA and r_quality[r] >= threshold
        if f_high and r_high and f_base == r_base:
            high_quality_matches += 1
            base = f_base
        elif (f_high and r_high) or (start <= column <= end and (f < 0 or r < 0)):
            base = "N"
        elif f_high:
            base = f_base
        elif r_high:
            base = r_base
        else:
            base = "N"
        consensus.append(base)

    problem = ""
    if overlap < min_overlap or identity < min_identity or high_quality_matches < min_high_quality_matches:
        problem = "Phần chồng lấn chưa đạt tiêu chí để ghép một consensus đáng tin cậy."
    else:
        alternative = next(alignments, None)
        if alternative is not None and alternative.coordinates.tolist() != alignment.coordinates.tolist():
            problem = "Có nhiều cách căn chỉnh đồng điểm. Chưa xác định được vị trí ghép duy nhất."
    return ReadAssembly(
        forward_aligned=f_aligned, reverse_aligned=r_aligned, overlap=overlap,
        identity=identity, high_quality_matches=high_quality_matches,
        sequence=None if problem else "".join(consensus), problem=problem,
    )
