# AFMのシミュレーション
## 共通
main.exe Probe_Model Surface_Model Scan_Type その他
### Surface_Model
- probe0 (単独原子)
- probe1 (半球)
- probe2 (直線 + 半球)
- probe3 (三角錐 + 半球)
- probe4 (CO + 半球)

### Probe_Model
- surface0 (単独原子)
- surface1 (2*2の原子)

### Scan_Type
- scan0 (xy平面)
- scan1 (z方向)


### 具体的な例
- 原子のポテンシャルの観測
main.exe probe0 surface0 scan1 0 0
- 曲率70+CO探針のxy(z=3Å)の観察画像
main.exe probe4 surface0 scan0 70 3
- 曲率70+高さ5層の直線(z=3Å)の観察画像
main.exe probe2 surface0 scan0 70 3 5

