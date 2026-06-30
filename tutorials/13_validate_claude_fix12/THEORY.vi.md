# LÝ THUYẾT — quá trình mô phỏng tổng hợp liên tục hoạt động ra sao, và nó mô hình hoá điều gì

> 🇬🇧 English version: [`THEORY.md`](THEORY.md).

Tài liệu này giải thích chi tiết phần vật lý và sinh học đằng sau lần chạy trong
tutorial này: cách `topo-csp` "mọc" dần một chuỗi protein non (nascent) trên ribosome,
và **ba tiểu giai đoạn (sub-stage) mô phỏng tương ứng với chu trình kéo dài (elongation
cycle) sinh học thật như thế nào**. Đây là tài liệu lý thuyết đi kèm với
[`README.md`](README.md) (phần thực hành), và **hợp nhất phần cơ chế từng giai đoạn**
vốn trước đây nằm trong một tệp `STAGES.md` riêng (nay đã gộp vào §3 bên dưới), nên tài
liệu này tự chứa (self-contained).

Tutorial 13 chạy trọn vẹn chuỗi **4c5c gồm 306 gốc axit amin (residue)** từ đầu đến
cuối, và là bước kiểm chứng (validation) rằng giao thức ổn định về mặt số học trên toàn
bộ một protein (xem §7).

---

## 1. Đối tượng được mô hình hoá: tổng hợp đồng dịch mã (co-translational synthesis)

Trong tế bào sống, protein không gập từ một chuỗi dài có sẵn. Nó được **tổng hợp theo
hướng, từ đầu N trước**, bởi ribosome, mỗi lần thêm một axit amin, trong khi chuỗi đang
lớn dần ("nascent") luồn ra ngoài qua **đường hầm thoát (exit tunnel)** của ribosome
(~80 Å chiều dài, ~10–20 Å chiều rộng) và bắt đầu **gập đồng dịch mã (co-translational
folding)** — tức là gập *ngay khi nó ló ra*. Vì vậy động học của quá trình tổng hợp rất
quan trọng: thời gian ribosome lưu lại (dwell) trên mỗi codon quyết định mỗi đoạn chuỗi
có bao nhiêu thời gian để lấy mẫu (sample) các cấu dạng trước khi gốc tiếp theo được thêm
vào. Các codon hiếm (được giải mã chậm) có thể đóng vai trò "điểm dừng dịch mã"
(translational pause) làm thay đổi kết quả gập.

Mô phỏng tái hiện quá trình này: nó cho một protein thô hạt (coarse-grained, CG) lớn dần
từng hạt một, đi ra khỏi một ribosome thô hạt, **định thời mỗi gốc theo codon của nó trên
mRNA**, để chuỗi non gập dưới động học thực tế, phân giải theo codon.

### Chu trình kéo dài thật (thêm một axit amin)

Quá trình kéo dài dịch mã ở vi khuẩn lặp lại một chu trình sinh hoá ba bước cho mỗi codon:

1. **Chọn aminoacyl-tRNA / giải mã (decoding).** Phức bậc ba (ternary complex,
   EF-Tu·GTP·aminoacyl-tRNA) đưa một aa-tRNA vào **vị trí A (A site)** của ribosome. Việc
   bắt cặp codon–anticodon đúng kích hoạt thuỷ phân GTP và "an vị" (accommodation) của
   aa-tRNA. Đây là bước **phụ thuộc codon, biến thiên mạnh, và thường là bước giới hạn
   tốc độ** — thời gian của nó phụ thuộc vào độ phong phú của tRNA tương ứng (codon usage
   bias).
2. **Chuyển peptidyl (tạo liên kết peptide).** Trung tâm peptidyl-transferase (PTC) chuyển
   chuỗi peptide non từ tRNA ở **vị trí P (P site)** sang aminoacyl-tRNA ở vị trí A. Chuỗi
   bây giờ dài thêm một gốc và gắn vào tRNA ở vị trí A. Nhanh (~0,3 ms).
3. **Dịch chuyển vị trí (translocation).** EF-G·GTP đẩy ribosome tiến tới một codon: các
   tRNA dịch **A→P** và **P→E**, tRNA đã khử acyl rời đi qua vị trí E, và vị trí A được
   giải phóng cho aa-tRNA kế tiếp. ~vài ms.

Giao thức Tổng hợp Liên tục (Continuous Synthesis Protocol, CSP;
`continuous_synthesis_v6.py`) của O'Brien chia thời gian lưu trên mỗi codon thành đúng ba
phần này, và `topo.csp` tái hiện chúng dưới dạng **ba tiểu giai đoạn MD cho mỗi gốc**.

---

## 2. Mô hình mô phỏng (một bước kéo dài)

Hệ mô phỏng được dựng từ bộ máy đã được tái sử dụng và kiểm thử
(`topo.csp.core` + `topo.csp.ribosome`):

- **Protein non — chuỗi thô hạt dựa-trên-cấu-trúc (structure-based / Gō-like).** Một hạt
  cho mỗi gốc, đặt tại vị trí Cα. Các tiếp xúc bản địa (native contacts — cặp gốc gần nhau
  trong cấu trúc tinh thể 4c5c đã gập) nhận giếng thế hút; mọi cặp khác là đẩy. **Do đó
  cấu trúc bản địa là cực tiểu năng lượng**, nên chuỗi đang lớn dần gập *về phía* cấu trúc
  bản địa của 4c5c — đồng dịch mã, đúng như chuỗi thật gập khi ló ra. Liên kết là **điều
  hoà mềm dẻo (flexible harmonic)**, không phải ràng buộc cứng — xem §7.
- **Ribosome — khung cảnh cứng (rigid scenery).** Ribosome 50S + tRNA đã cắt cụt, dạng CG
  (`ribosome_trunc.pdb`, ~4.600 hạt khối lượng 0) được cố định trong không gian. Nó không
  chuyển động, nhưng **tương tác thể tích loại trừ (excluded volume) và tĩnh điện với
  chuỗi non vẫn được bật**, nên chuỗi "cảm nhận" được thành đường hầm và bề mặt ribosome.
  Trục đường hầm trùng với **+x** (chuỗi thoát ra theo hướng +x).
- **Các điểm neo tại PTC.** Hai hạt ribosome được chọn làm điểm tham chiếu cố định: **điểm
  neo P (P-anchor)** (hạt "R" của gốc 76 của tRNA ở vị trí P) và **điểm neo A (A-anchor)**
  (hạt "R" của gốc 76 của tRNA ở vị trí A). Chúng đại diện cho nơi mà peptidyl-tRNA (P) và
  aminoacyl-tRNA đi vào (A) giữ đầu C của chuỗi.
- **Dây neo đầu C — một ràng buộc điều hoà (harmonic restraint).** Hạt đầu C hiện thời bị
  ràng buộc về một trong các điểm neo với `U = k·|r − r₀|²`, `k = restraint_k = 83680
  kJ/mol/nm²` (= 200 kcal/mol/Å²). Điều này tái hiện sự gắn cộng hoá trị của đầu C chuỗi
  non vào tRNA đang nằm ở vị trí A hoặc P. **Việc chuyển mục tiêu ràng buộc từ điểm neo A
  sang điểm neo P chính là cách tái hiện quá trình translocation.**
- **Thành đường hầm — một mặt phẳng một phía (one-sided plane).** Một thành nửa-điều-hoà
  giữ mọi hạt non ở `x ≥ x₀` (`k = 8368 kJ/mol/nm²`). Vì 50S đã bị **cắt cụt** cho nhanh,
  không có hạt ribosome nào bên dưới PTC — phần hầm dưới là khoảng trống; thành tường cung
  cấp "sàn" còn thiếu để chuỗi chỉ có thể đùn ra **phía trước** (+x) và không tụt ngược
  qua điểm tổng hợp vào vùng đã cắt cụt. Mặt phẳng `x₀` được **tự động suy ra từ cấu trúc
  ribosome** (không phải tham số cấu hình): mặt giữ đầu C ở vị trí P (thấp hơn),
  `x₀ = min(P-anchor.x, A-anchor.x) + ptc_offset` — với `ribosome_trunc.pdb` này là
  `0,5705 + 0,476 = 1,0465 nm` (≈ giá trị cũ 1,05 đặt cứng).
- **Bộ điều nhiệt (thermostat).** Động lực học Langevin tại `ref_t = 310 K`, ma sát
  `tau_t = 0,01 /ps`, bước thời gian `dt = 0,015 ps`.

### Tính tiếp xúc một-lần-rồi-cắt (build-once-subset)

Bản đồ tiếp xúc được tính **một lần** trên toàn cấu trúc 4c5c (`R_full`, `eps_full`,
306×306). Với chiều dài non `L`, mô hình dùng khối con `L×L` ở góc trên-trái — STRIDE và
phép phân tích tiếp xúc theo nguyên tử nặng không bao giờ bị chạy lại cho mỗi chiều dài.
Nhờ vậy, các gốc `1..L` mang đúng những tiếp xúc bản địa mà chúng sẽ có trong cấu trúc gập
cuối cùng, và chuỗi có thể hình thành cấu trúc bản địa ngay khi các gốc liên quan tồn tại.

---

## 3. Ba giai đoạn: sinh học ↔ mô phỏng

Mỗi axit amin được thêm qua **một chu trình kéo dài ribosome**, được O'Brien chia thành
ba tiểu bước động học; `topo-csp` chạy **một đoạn MD cho mỗi tiểu bước**. Với chiều dài non
`L`, mỗi tiểu giai đoạn là một mô phỏng ngắn độc lập (thư mục riêng
`synth_out/L_<L>/stage_<s>/`); cấu trúc cuối của giai đoạn 3 làm hạt giống (seed) cho giai
đoạn 1 của gốc kế tiếp. Tổng thời gian dịch mã in-vivo của codon được **chia cho ba giai
đoạn**, và vòng lặp gói gọn trong một dòng:

```
codon_total  =  GĐ 1       +   GĐ 2               +   GĐ 3
                (cố định)      (cố định + traffic)    (phần dư, phụ thuộc codon)

đưa vào A (GĐ 1) → giữ ở A trong khi translocation (GĐ 2)
   → chuyển sang P và chờ tRNA kế tiếp (GĐ 3) → lặp lại cho gốc L+1
```

| GĐ | Quá trình sinh học thật | Mô phỏng làm gì | Đầu C ràng buộc về | Thời gian lưu trung bình |
|----|--------------------------|------------------|--------------------|---------------------------|
| **1** | **Chuyển peptidyl** — tạo liên kết peptide; chuỗi nay nằm trên tRNA ở vị trí A | Hạt mới `L` được **đặt tại điểm neo A** (+`buffer = 0,4 nm` vào trong hầm), liên kết với gốc `L−1`; tối thiểu hoá năng lượng; chạy MD | **điểm neo A** | `time_stage_1 = 0,34 ms` |
| **2** | **Bắt đầu translocation** — EF-G bắt đầu đẩy ribosome tiến tới | Tiếp tục từ GĐ1, **vẫn giữ ở điểm neo A**; chạy MD | **điểm neo A** | `time_stage_2 = 4,20 ms` |
| **3** | **Hoàn tất translocation + chờ aa-tRNA kế tiếp** — peptidyl-tRNA nay ở vị trí P; ribosome chờ giải mã codon kế tiếp | **Chuyển ràng buộc A→P** (chính cú dịch chuyển hình học này *là* translocation), rồi chạy MD để chuỗi giãn/gập | **điểm neo P** | phần dư = (tổng thời gian codon kế tiếp) − GĐ1 − GĐ2 |

### Giai đoạn 1 — chuyển peptidyl / đưa vào vị trí A
- **Sinh học.** tRNA của axit amin đi vào được an vị tại vị trí A và PTC tạo liên kết
  peptide mới. Chuỗi non nay dài `L`, được giữ cộng hoá trị bởi tRNA ở vị trí A.
- **Mô phỏng.** Hạt `L` được gieo tại điểm neo A (lệch `buffer` theo +x để không chồng lên
  hạt neo). Một liên kết mềm dẻo tới gốc `L−1` được tạo ra — ban đầu nó bị kéo giãn (điểm
  neo A và P cách nhau ~0,92 nm), nên cấu trúc được **tối thiểu hoá năng lượng** trước để
  giãn vị trí, sau đó MD chạy theo số bước của GĐ1. Đầu C bị ràng buộc về điểm neo A.
- **Thời gian trung bình.** `time_stage_1 = 0,000340 s` (thời gian lưu thực nghiệm của việc
  tạo liên kết peptide).

### Giai đoạn 2 — bắt đầu translocation
- **Sinh học.** EF-G·GTP gắn vào và ribosome bắt đầu di chuyển một codon dọc mRNA; các tRNA
  bắt đầu dịch A→P / P→E.
- **Mô phỏng.** Hệ tiếp tục từ cấu trúc cuối của GĐ1, **vẫn ràng buộc tại điểm neo A**, và
  chạy MD theo số bước của GĐ2. (Cú dịch chuyển hình học A→P thực sự được áp dụng ở đầu
  GĐ3.)
- **Thời gian trung bình.** `time_stage_2 = 0,004201 s` (thời gian lưu translocation thực
  nghiệm), cộng thêm hiệu chỉnh **giao thông ribosome (ribosome traffic)** tuỳ chọn (tắt
  trong tutorial này; xem §5).

### Giai đoạn 3 — dịch chuyển sang vị trí P + chờ gắn tRNA
- **Sinh học.** Translocation hoàn tất (peptidyl-tRNA nay ở **vị trí P**, vị trí A trống) và
  ribosome **chờ aminoacyl-tRNA tương ứng kế tiếp** đến và được giải mã. Sự chờ đợi này
  **phụ thuộc codon** và thường là phần dài nhất của chu trình — chính nó làm các codon
  hiếm trở nên chậm.
- **Mô phỏng.** Ràng buộc đầu C **chuyển từ điểm neo A sang điểm neo P** — tái hiện cú
  translocation một-codon A→P dưới dạng một cú dịch chuyển hình học rời rạc — và MD chạy
  theo số bước của GĐ3. Với đầu C bị kéo về vị trí P và thành đường hầm ngăn gập ngược,
  chuỗi đùn thêm một nấc nữa về phía lối ra. Cấu trúc cuối của GĐ3 trở thành **hạt giống cho
  GĐ1 của gốc `L+1`**.
- **Thời gian trung bình.** phần **dư**: (thời gian in-vivo trung bình của codon *kế tiếp*)
  − `time_stage_1` − `time_stage_2`. Nếu một codon nhanh khiến giá trị này ≤ 0 thì nó được
  đặt sàn về một số dương rất nhỏ. Đây là nơi tính biến thiên theo từng codon đi vào lịch
  trình.

> **Cơ chế và định thời — một lưu ý trung thực.** Việc chuyển ràng buộc (cú dịch chuyển
> hình học A→P tức thời) xảy ra ở **đầu GĐ3**, trong khi **thời lượng** được gán cho
> translocation là **GĐ2** và thời lượng gán cho việc chờ giải mã là **GĐ3**. Do đó cú dịch
> chuyển vật lý và khoảng thời gian mang nhãn "translocation" bị tách rời nhẹ — đây là một
> sự đơn giản hoá có chủ đích. Tương tự, liên kết peptide hiện diện trong mô hình liên kết
> ngay từ GĐ1 thay vì được "bật" giữa giai đoạn, và hình học liên kết của tRNA tại A/P
> không được mô hình hoá tường minh. Phần **định thời** (ba thời gian lưu phân giải theo
> codon cho mỗi gốc) trung thành với O'Brien; phần **cơ chế** của từng giai đoạn là một mô
> hình rút gọn (hình học liên kết tRNA tại A/P không được mô hình hoá tường minh).

---

## 4. Từ codon đến số bước MD (động học)

Lõi định thời là `topo.csp.kinetics` (thuần tuý, không có OpenMM). Với mỗi gốc nó trả lời:
*mỗi tiểu giai đoạn chạy bao nhiêu bước tích phân?*

**(a) Thời gian dịch mã trung bình theo codon.** mRNA (`4c5c_mrna.txt`) được tách thành các
codon; một bảng tra (`trans_times.txt`, bảng Fluitt *E. coli* tại 310 K) ánh xạ mỗi codon
sang **thời gian dịch mã in-vivo trung bình** của nó tính bằng giây. Đây là **thời gian
chuyển tiếp lần đầu trung bình (mean first-passage time, mFPT)** nội tại của codon — chi
phối bởi độ sẵn có của tRNA tương ứng (codon usage). Gọi là `τ(codon)`.

**(b) Lấy mẫu thời gian chuyển tiếp lần đầu (FPT).** Quá trình kéo dài thật mang tính ngẫu
nhiên: mỗi bước bị chi phối bởi một sự kiện phân tử giới hạn tốc độ, nên thời gian chờ của
nó tuân theo **phân phối mũ (exponential)**. Thời gian lưu thực tế của mỗi tiểu giai đoạn
được lấy mẫu theo `t = −mean · ln(U)`, `U ∼ Uniform(0,1)` (`random.expovariate`). Một
`random_seed` cố định khiến lịch trình có thể tái lập.

**(c) Phép chia ba giai đoạn** cho chiều dài non `L` (đánh chỉ số từ 1):
```
t1  ~ Exp(  time_stage_1 )                              # chuyển peptidyl (trung bình cố định)
t2  ~ Exp(  time_stage_2 + max(0, hiệu_chỉnh_traffic) ) # translocation (+ traffic tuỳ chọn)
t3  ~ Exp(  τ(codon kế tiếp) − time_stage_1 − time_stage_2 )  # chờ giải mã codon kế tiếp
```
nên **đồng hồ mỗi chu trình** = liên kết peptide cố định + translocation cố định + chờ giải
mã biến thiên. (Lưu ý chỉ số: GĐ3 dùng thời gian trung bình của codon *kế tiếp* — ribosome,
sau khi vừa lắp gốc `L`, nay chờ tRNA cho gốc `L+1`.)

**(d) Giây in-vivo → số bước in-silico.** Mô hình thô hạt tiến triển trên thang thời gian
ngắn hơn rất nhiều so với dịch mã thật, nên một **hệ số nén thời gian (time-compression
factor)** ánh xạ giây sang số bước MD:
```
t_sim (ns)  =  t_s · 1e9 / scale_factor
n_steps     =  t_sim (ns) / dt(ns) ,   dt(ns) = dt_ps · 1e-3
```
**`scale_factor` lớn hơn ⇒ ít bước hơn cho mỗi gốc ⇒ chạy nhanh hơn**, trong khi vẫn giữ
được định thời *tương đối* giữa codon nhanh và codon chậm. Tutorial này dùng
`scale_factor = 216564650` (= **50×** giá trị sản xuất 4331293) để toàn bộ chuỗi 306 gốc
hoàn tất nhanh; *tỷ lệ* giữa các codon — phần vật lý quan trọng cho gập đồng dịch mã — không
đổi. Số bước còn bị kẹp (clamp) về `[min_steps_per_stage, max_steps_per_stage] = [400,
10000]` cho khả thi (chỉ kẹp số bước MD; **thời gian lưu (giây)** đã lấy mẫu vẫn được ghi
nguyên vẹn trong `synth_out/dwell_times.dat`).

---

## 5. Giao thông ribosome, và hai loại thời gian theo codon: `intrinsic` và `real`

Phần động học mang **hai** danh sách thời gian theo codon:

- **`intrinsic[i]`** — thời gian dịch mã trung bình "trần" của codon, tra thẳng từ bảng
  `trans_times` (không xếp hàng). Đây chính là mFPT nội tại ở §4.
- **`real[i]`** — cũng thời gian đó nhưng **đã hiệu chỉnh giao thông ribosome (traffic)**.
  Trong tế bào thật, các ribosome xếp hàng trên cùng một mRNA: một ribosome chậm phía trước
  làm trễ những ribosome phía sau, kéo dài thời gian lưu hiệu dụng trên mỗi codon. O'Brien mô
  hình hoá điều này bằng chương trình ngoài `ribosome_traffic` — cho biết mRNA, các thời gian
  nội tại và `initiation_rate`, nó trả về thời gian theo codon đã hiệu chỉnh traffic.

Chỉ **GĐ2** dùng phần hiệu chỉnh này:
```
correction        = real[L−1] − intrinsic[L−1]
mean(GĐ 2)        = time_stage_2 + correction   nếu correction > 0,   ngược lại = time_stage_2
```

**Khi `ribosome_traffic = no`, `real == intrinsic`.** `build_mfpt_lists` khởi tạo
`real = list(intrinsic)` và **chỉ** ghi đè khi traffic vừa được yêu cầu *vừa* chạy được
chương trình ngoài; nếu không, `real` được giữ bằng `intrinsic`. Vậy `real == intrinsic` mỗi
khi:

- `ribosome_traffic = no` (nhánh `if ribosome_traffic:` bị bỏ qua hoàn toàn), **hoặc**
- `ribosome_traffic = yes` nhưng chương trình `ribosome_traffic` thiếu / lỗi (nó trả về
  `None` và mã quay về `real = intrinsic`, kèm một cảnh báo in ra).

Trong mọi trường hợp đó `correction = 0`, nên **trung bình của GĐ2 đơn giản bằng
`time_stage_2`** (không bị kéo dài). Tutorial này đặt `ribosome_traffic = no`, nên
`real == intrinsic` ở mọi nơi và lịch trình được quyết định hoàn toàn bởi các thời gian theo
codon "trần"; ảnh hưởng của traffic ở các chiều dài này là nhỏ.

Mã nguồn: `topo.csp.kinetics.build_mfpt_lists` / `ribosome_traffic_times` / `stage_dwell_times`.

---

## 6. Sau gốc cuối cùng: phóng thích (ejection) (và phân ly)

Khi gốc 306 được thêm vào, protein đã hoàn chỉnh. Mô phỏng sau đó chạy một **giai đoạn
phóng thích sau tổng hợp** (`ejection_steps = 50000`): ràng buộc đầu C được **giải phóng**
(cắt dây neo), trong khi ribosome cứng và thành đường hầm một phía vẫn còn. Về mặt sinh học
đây là **kết thúc (termination)** — các yếu tố giải phóng (release factors) thuỷ phân liên
kết peptidyl-tRNA và protein hoàn chỉnh được phóng thích. Khi dây neo biến mất, chuỗi
**khuếch tán ra khỏi đường hầm dọc theo +x** (thành một phía đóng vai trò rào phản xạ làm
lệch chuyển động về phía trước) và rời khỏi ribosome, tự do hoàn tất việc gập. Một giai đoạn
**phân ly (dissociation)** tuỳ chọn (`dissociation_steps = 0` ở đây) sẽ tiếp tục cho protein
tự do đi xa khỏi ribosome.

Để có một minh hoạ thoát ra dài hơi, chuyên biệt, xem `eject_demo.py` (giải phóng chuỗi cuối
và chạy đủ lâu để nó rời hẳn khỏi đường hầm).

---

## 7. Tích phân số học và "lớp bảo vệ ổn định" (vì sao Tutorial 13 tồn tại)

Chuỗi được tích phân với **liên kết điều hoà mềm dẻo** tại `dt = 0,015 ps`. Cần liên kết
mềm dẻo (thay vì ràng buộc cứng `AllBonds`) vì GĐ1 gieo hạt mới cách bạn liên kết của nó
~1 nm (đưa vào vị trí A), điều mà một ràng buộc khoảng cách cứng không thể biểu diễn — một
liên kết điều hoà hấp thụ độ giãn và bộ tối thiểu hoá làm giãn nó ra.

Cái giá phải trả: tại 15 fs, việc tích phân chỉ **ổn định ở mức biên** đối với một số cấu
hình. Khi một gốc vừa được thêm tạo nên một **tiếp xúc bản địa (Gō) cứng**, chu kỳ dao động
của tiếp xúc đó giảm xuống dưới mức mà một bước 15 fs có thể tích phân, và động lực học
**phân kỳ** (thế năng → ~10¹³ kJ/mol), làm hỏng các khung hình (frame) của giai đoạn đó. Đây
chính là lỗi tiềm ẩn trong các lần chạy đầy đủ chưa có lớp bảo vệ: ~5 trên 306 giai đoạn nổ tung
trong một lần chạy đầy đủ. Lỗi này **mang tính tất định theo bước thời gian, không ngẫu
nhiên** — đo trên 4c5c, chiều dài 10 phân kỳ ở 15 fs với mọi hạt giống vận tốc nhưng ổn định
ở 7,5 fs với mọi hạt giống. (Bản tham chiếu của O'Brien tránh hẳn vấn đề bằng cách dùng ràng
buộc cứng `AllBonds`, loại bỏ mode liên kết nhanh.)

**Cách khắc phục** (trong `topo.csp.core.run_length`, lớp bảo vệ ổn định mỗi
giai đoạn): mỗi giai đoạn được chạy theo từng khúc (chunk) trong khi theo dõi **giá trị cực
đại** của |PotE|; nếu một giai đoạn phân kỳ (max |PotE| > 10⁹ kJ/mol) nó được **chạy lại với
bước thời gian giảm một nửa và số bước nhân đôi**. Vì thời gian lưu vật lý là `n_steps · dt`,
việc giảm `dt` một nửa và nhân đôi `n_steps` **giữ nguyên chính xác thời gian lưu** trong khi
làm cho việc tích phân ổn định (tối đa 6 lần giảm nửa). Trường hợp thông thường chỉ chạy một
lần ở 15 fs.

**Tutorial 13 chính là phép kiểm chứng cách khắc phục đó ở quy mô đầy đủ.** Khi tổng hợp cả
306 gốc (919 giai đoạn), **không** giai đoạn nào nổ tung và max |PotE| tệ nhất trên toàn bộ
lần chạy chỉ ~1,8×10³ kJ/mol — xác nhận rằng lớp bảo vệ duy trì được trên toàn chuỗi, đúng
trong chính chế độ mà lần chạy đầy đủ chưa có lớp bảo vệ đã thất bại.

---

## 8. Hệ 4c5c và các tham số của tutorial này

- **Protein:** 4c5c, **306 gốc**, ba miền (domain) (từ `domain.yaml`, ánh xạ từ mô hình CG
  của O'Brien): **A** = 1–84, **B** = 85–110 + 184–306 (không liền mạch), **C** = 111–183,
  với hệ số nhân độ mạnh tiếp xúc Gō theo từng miền và từng giao diện (bản tương tự
  dựa-trên-cấu-trúc của `nscal` trong O'Brien).
- **mRNA:** `4c5c_mrna.txt` — mỗi gốc một codon (+ codon dừng), **giống hệt từng byte** với
  mRNA "fast" tham chiếu của O'Brien, nên lịch trình codon đúng bằng của bản tham chiếu.
- **Các con số then chốt** (`csp_val.ini`): `dt = 0,015 ps`, `ref_t = 310 K`,
  `scale_factor = 216564650` (50× sản xuất), `time_stage_1 = 0,34 ms`,
  `time_stage_2 = 4,20 ms`, số bước kẹp về `[400, 10000]`, `restraint_k = 83680`, thành hầm
  tại `x ≥ 1,05 nm` (`k = 8368`), `ejection_steps = 50000`.

> **Lưu ý về tệp cấu hình.** Phần chú thích ở đầu `csp_val.ini` đã cũ (sao chép từ Tutorial
> 12 — nhắc đến "length 1 → 10" và `scale_factor = 4331293`). Các thiết lập *đang hoạt động*
> là các giá trị nêu trên (`L_max = 306`, `scale_factor = 216564650`); hãy tin các khoá
> (key), đừng tin khối chú thích.

---

## 9. Tóm tắt: bảng tương ứng nhanh

```
CHU TRÌNH KÉO DÀI THẬT (mỗi codon)          MÔ PHỎNG (mỗi gốc L, 3 tiểu giai đoạn MD)
─────────────────────────────────          ───────────────────────────────────────────
tạo liên kết peptide (PTC)          →  GĐ 1: đặt hạt L tại điểm neo A, giữ ở A, chạy MD
                                              trung bình = time_stage_1 (0,34 ms)
EF-G bắt đầu translocation           →  GĐ 2: tiếp tục, vẫn giữ ở A, chạy MD
                                              trung bình = time_stage_2 (4,20 ms)
translocation hoàn tất (→ vị trí P)  →  GĐ 3: chuyển ràng buộc A→P (cú dịch chuyển), MD
+ chờ giải mã aa-tRNA kế tiếp                 trung bình = thời_gian_codon_kế − GĐ1 − GĐ2
(phụ thuộc codon, giới hạn tốc độ)            (đây là phần biến thiên, phân giải theo codon)

kết thúc / phóng thích               →  ejection: cắt dây neo, chuỗi khuếch tán ra +x
```

Các thời gian lưu phân giải theo codon (lấy mẫu mũ, nén bởi `scale_factor`) quyết định mỗi
đoạn chuỗi giãn trong bao lâu trước khi gốc kế tiếp đến — nhờ vậy chuỗi gập **dưới cùng động
học tương đối như dịch mã thật**, trong khi khung cảnh ribosome + thành đường hầm giữ cho nó
luồn về phía trước ra khỏi đường hầm thoát.

---

## Tham chiếu mã nguồn & tài liệu

- Động học: [`topo/csp/kinetics.py`](../../topo/csp/kinetics.py) — `codon_mfpt_list`,
  `sample_fpt`, `stage_dwell_times`, `seconds_to_steps`, `stage_steps`.
- Vòng lặp 3 giai đoạn: [`topo/csp/protocol.py`](../../topo/csp/protocol.py) —
  `run_continuous_synthesis` (ba lần gọi `run_length` cho mỗi gốc; chuyển ràng buộc A→A→P;
  giai đoạn ejection).
- MD mỗi giai đoạn + lớp bảo vệ ổn định:
  [`topo/csp/elongate.py`](../../topo/csp/elongate.py) — `run_length`.
- Cơ chế ba giai đoạn nằm ở §3 của tài liệu này (hợp nhất từ tệp `STAGES.md` cũ); cách
  khắc phục số học và lý do ở §7; các con số kiểm chứng ở [`README.md`](README.md).
