# Hong 2017 - J. Phys. D

## 论文信息
- **标题**: Plasma-catalytic synthesis of ammonia in a dielectric barrier discharge
- **期刊**: Journal of Physics D: Applied Physics, 50, 154005 (2017)
- **机制**: 基线 N2-H2 机制，14个反应类别

## kinet.inp 状态
- ✅ **已通过 ZDPlasKin preprocessor 验证**
- 元素：E, N, H, S, M (5个)
- 物种：35个
- 反应：56个（21个 BOLSIG+ 电子碰撞反应）
- 文件大小：~89行

## 已知问题
- 物种 N2H, N(2D), N(2P), H^- 声明但未使用
- 需要 BOLSIG+ 截面数据（`bolsigdb.dat`）才能运行模拟
- 缺少完整的 Tables 1-3 数据（只有4页 ISPC 预印本）

## 下一步
1. 获取 N2/H2 截面数据
2. 运行 `main_bolsig.exe` 验证完整模拟
3. 获取完整论文以核对反应速率系数
