# 主要文件
1. copper{923,973,1023}.i
2. exp{923,973,1023}K.csv, time(s)--gr0_area($\mu m^2$)
3. simu_copper{923,973,1023}K.csv, time(s)--gr0_area($\times 10^{-2} \mu m^2$)

# 模型参数

## 第一级参数
1. GB mobility perfactor, $m^4/(J \cdot s)$
2. GB energy, $J/m^2$, $0.0708$
3. GB width, $\times 10^{-7} m$, $1$

## 二级参数：相场模型模型参数
$$
\begin{aligned}
	\gamma &=1.5\\
	\kappa _i&\approx \frac{3}{4}\sigma _{\mathrm{GBB}}\omega _{\mathrm{GB}}\\
	L&\approx \frac{4}{3}\frac{m_{\mathrm{GB}}}{\omega _{\mathrm{GB}}}\\
	\mu &\approx \frac{3}{4}\frac{1}{f_{0, \mathrm{sadde}}(\gamma )}\frac{\sigma _{\mathrm{GB}}}{\omega _{\mathrm{GB}}}\\
\end{aligned}
$$

# 模型参数优化过程
~~