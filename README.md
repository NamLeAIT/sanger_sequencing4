# Phân tích Sanger

App Streamlit cho dữ liệu Sanger dạng AB1. Tải lên một hoặc hai tệp Forward và Reverse, nhập tham chiếu nếu cần rồi bấm **Phân tích**.

- **Tải FASTA F/R** có nút tải riêng cho từng chiều ngay sau khi phân tích. Chọn một read để xem bảng từng base và tải CSV của read đó.
- **Chất lượng** hiển thị số base gọi được, tỷ lệ đạt ngưỡng Q, trung vị và đoạn đạt ngưỡng liên tục dài nhất. Ngưỡng Q nằm trên các tab kết quả.
- **Chromatogram** hiển thị bốn kênh huỳnh quang đã xử lý. Rê chuột lên điểm gọi base để xem vị trí và điểm Q.
- **So sánh tham chiếu** nhận chuỗi DNA hoặc FASTA. App chọn hướng có điểm căn chỉnh cục bộ cao hơn cho từng read. Kết quả hiển thị từng base, gap và tọa độ trong tệp gốc.
- **Consensus F/R** xuất hiện khi tải đủ F và R và để trống ô tham chiếu. App tìm phần chồng lấn giữa hai read và xuất FASTA nếu phần chồng lấn đạt tiêu chí.
- Trong mục **Kiểm tra từng vị trí**, chọn một tọa độ tham chiếu để xem base call, Q và peak tương ứng. Tải CSV hoặc FASTA từ phần kết quả.

## Chạy trên máy tính

Cần Python 3.12. Mở terminal trong thư mục đã giải nén rồi chạy:

```sh
python -m venv .venv
```

PowerShell trên Windows:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

macOS hoặc Linux:

```sh
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run app.py
```

## Triển khai trên Streamlit Community Cloud

1. Đưa nội dung thư mục đã giải nén lên một repository GitHub. Đặt `app.py`, `sanger.py`, `requirements.txt` và thư mục `.streamlit` ở thư mục gốc của repository.
2. Mở [Streamlit Community Cloud](https://share.streamlit.io/) và chọn **Create app**. Chọn repository, nhánh và **Main file path** là `app.py`.
3. Chọn Python **3.12** trong **Advanced settings** nếu cần, rồi bấm **Deploy**.

Gói mã nguồn chỉ chứa app. Bạn có thể tải tệp AB1 lên qua giao diện sau khi app chạy. App xử lý đầu vào trong phiên Streamlit và không ghi dữ liệu mẫu vào repository.

## Diễn giải kết quả

FASTA của một read giữ thứ tự và toàn bộ tọa độ máy đọc. Base A/C/G/T có điểm Q đạt ngưỡng được giữ lại. Vị trí điểm Q thấp, thiếu Q hoặc base mơ hồ nhận N. Có thể xuất từng read độc lập khi dùng một tệp AB1.

App dùng căn chỉnh cục bộ có gap. Tỷ lệ identity tính trên các cột của đoạn căn chỉnh, còn độ phủ tính theo toàn bộ độ dài tham chiếu. Nếu hai hướng đồng điểm, app ưu tiên chiều khai báo F/R. Các ký hiệu mơ hồ như N được hỗ trợ trong tham chiếu và được đánh dấu bằng dấu hỏi trong căn chỉnh.

Consensus dùng tọa độ của tham chiếu và chỉ lấy base có điểm Q đạt ngưỡng. Khi hai chiều có base đạt ngưỡng nhưng bất đồng hoặc không có base nào đạt ngưỡng, vị trí đó nhận N. FASTA consensus biểu diễn các vị trí trên tham chiếu. Những base bị căn chỉnh thành chèn xuất hiện thành dòng riêng trong CSV.

Khi không có tham chiếu, app cắt các đầu read không đạt ngưỡng Q, lấy đảo bổ sung của R và căn chỉnh F với R. Mặc định cần ít nhất 20 cặp base chồng lấn, 85% cột khớp trong đoạn chồng lấn và 8 vị trí mà cả hai read đạt ngưỡng Q và gọi cùng base. Bạn có thể chỉnh ba tiêu chí trong tab **Consensus F/R**. Nếu có nhiều căn chỉnh tối ưu hoặc phần chồng lấn dưới tiêu chí, app hiển thị lý do để bạn kiểm tra hai read. FASTA từ hai read giữ phần kéo dài đạt Q và gán N cho vị trí thiếu bằng chứng, xung đột chất lượng cao hoặc có gap ở giữa đoạn chồng lấn.
