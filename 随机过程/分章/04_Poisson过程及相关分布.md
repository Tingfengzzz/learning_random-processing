# 第04章 齐次Poisson过程及相关分布

来源：C3 p1–38；对应视频5。课件要求：掌握Poisson过程的定义及数字特征；理解与指数、Gamma、均匀分布的联系；熟练使用发生时刻、时间间隔的结论。

## 4.1 从Poisson分布到计数过程

【课件 C3 p2–4】若 $X\sim\operatorname{Poi}(a)$，

$$P(X=k)=e^{-a}\frac{a^k}{k!},\quad k=0,1,\ldots,\qquad E(X)=D(X)=a.$$

课件回顾的Poisson极限定理：当$n\to\infty$、$p=\lambda/n$时，固定非负整数$k$有

$$\binom nkp^k(1-p)^{n-k}\longrightarrow e^{-\lambda}\lambda^k/k!.$$

分布描述一个计数；过程则描述随时间变化的一族计数。计数过程 $N(t)$ 记录从0到$t$某类事件发生的次数，取非负整数值、不减，且 $N(t)-N(s)$ 记录 $(s,t]$ 的次数。

## 4.2 定义一：独立增量与Poisson分布增量

【课件 C3 p4】强度/率为 $\lambda>0$ 的齐次Poisson过程满足：

1. $N(0)=0$。
2. 不重叠区间的增量相互独立。
3. 对 $s\ge0,h>0$，$N(s+h)-N(s)\sim\operatorname{Poi}(\lambda h)$。

因此

$$P\{N(t)-N(s)=k\}=e^{-\lambda(t-s)}\frac{[\lambda(t-s)]^k}{k!},\quad0\le s<t.$$

**解题第一步永远是区间长度：** 率$\lambda$的单位若为“每小时”，时间就必须换成小时。Poisson分布参数是 $\lambda\times\text{时长}$，不是总把参数写成$\lambda$。

例如C3 p34流星例：平均每小时3颗，8点到12点长度4小时，无流星概率为 $e^{-12}$，不是$e^{-3}$。

## 4.3 定义二：平稳独立增量与小区间条件

【课件 C3 p6】除了 $N(0)=0$、平稳独立增量，还满足当 $h\downarrow0$ 时

$$P\{N(t+h)-N(t)\ge1\}=\lambda h+o(h),$$
$$P\{N(t+h)-N(t)\ge2\}=o(h).$$

这里 $o(h)/h\to0$。“发生至少一次”的主项与时长成正比；同一很短区间出现两次或更多次的概率比$h$更小。两式结合可得：恰好一次概率为 $\lambda h+o(h)$，零次概率为 $1-\lambda h+o(h)$。

课件用平稳性、稀有性、相继性解释这些条件。会解释含义、会与定义一对应即可，不擅自增加微分方程求解要求。

## 4.4 数字特征与相关性

【课件 C3 p5】由 $N(t)\sim\operatorname{Poi}(\lambda t)$，

$$E[N(t)]=D[N(t)]=\lambda t.$$

若$s\le t$，写 $N(t)=N(s)+[N(t)-N(s)]$；两项独立，因此

$$E[N(s)N(t)]=E[N(s)^2]+E[N(s)]E[N(t)-N(s)]$$
$$=\lambda s+(\lambda s)^2+\lambda s\lambda(t-s)=\lambda s+\lambda^2st.$$

对任意$s,t$，

$$r_N(s,t)=\lambda^2st+\lambda\min(s,t),\qquad R_N(s,t)=\lambda\min(s,t).$$

特征函数：$\phi_{N(t)}(u)=\exp\{\lambda t(e^{iu}-1)\}$。

**易错点：** $N(s)$和$N(t)$一般不独立，因为包含同一段历史计数；独立的是 $N(s)$ 与 $N(t)-N(s)$。Poisson过程具有平稳增量，但均值$\lambda t$变化，过程本身不是宽平稳。

## 4.5 分解/稀疏化与叠加

【课件 C3 p7–14】

### 独立分流

原过程率为$\lambda$。每个到达独立地以概率$p$进入系统Y、以$1-p$进入Z，分流选择与原过程独立。则Y的过程率为$\lambda p$，Z的过程率为$\lambda(1-p)$。

固定$t$，给定原计数 $N(t)=m$，有 $Y(t)\mid N(t)=m\sim\operatorname{Bin}(m,p)$。全概率展开：

$$P(Y(t)=k)=\sum_{m=k}^{\infty}\binom mkp^k(1-p)^{m-k}e^{-\lambda t}\frac{(\lambda t)^m}{m!}$$
$$=e^{-\lambda pt}\frac{(\lambda pt)^k}{k!}.$$

课件推广到$n$类，固定分流概率$p_k$满足和为1，各子过程率为$\lambda p_k$。在独立分配的标准条件下，分出的过程也相互独立；仅知道各自边缘是Poisson不足以断言这一点。

若分流概率随到达时刻变化为$p_k(s)$，正确的累计均值是

$$m_k(t)=\lambda\int_0^tp_k(s)\,ds.$$

【核对修正】C3 p10用平均概率 $\bar p_k(t)=t^{-1}\int_0^tp_k(s)ds$ 表示固定时刻的均值 $\lambda t\bar p_k(t)$。不能把依赖$t$的平均值误当成全程恒定强度；通常对应第05章的非齐次过程。

**课件应用：** 星球数量$\operatorname{Poi}(\lambda V)$、每个独立有生命概率$p$，有生命星球数为 $\operatorname{Poi}(\lambda pV)$（p11选A）；蚕产卵过程率$\lambda$、每个卵独立成活概率$p$，$[0,t]$成活$k$只概率为 $e^{-\lambda pt}(\lambda pt)^k/k!$。

### 独立叠加

若若干Poisson过程相互独立，率分别$\lambda_1,\ldots,\lambda_m$，总过程为Poisson过程，率是 $\sum\lambda_i$。

原课件汽车例：三种颜色汽车到达率$\lambda_1,\lambda_2,\lambda_3$，合并后计数服从 $\operatorname{Poi}[(\lambda_1+\lambda_2+\lambda_3)t]$。

相减一般不再是Poisson过程。独立$X,Y$的差 $Z=X-Y$ 可能为负，且不单调；其均值 $(\lambda_1-\lambda_2)t$ 与方差 $(\lambda_1+\lambda_2)t$ 通常不同。课件只需据此否定，不增加该差分布的名称和理论。

## 4.6 指数分布和Gamma分布回顾

【课件 C3 p16–23】指数分布采用率参数：

$$f_X(x)=\lambda e^{-\lambda x}\ (x>0),\quad F_X(x)=1-e^{-\lambda x}\ (x\ge0).$$

均值$1/\lambda$、方差$1/\lambda^2$。**无记忆性：**

$$P(X>s+t\mid X>s)=\frac{e^{-\lambda(s+t)}}{e^{-\lambda s}}=e^{-\lambda t}=P(X>t).$$

意思是已等$s$且事件未发生时，再等超过$t$的概率仍是原来的尾概率；不是说总寿命与已用时间相互独立。

Gamma分布（形状$a>0$、率$\lambda>0$）的密度为

$$f_X(x)=\frac{\lambda^a}{\Gamma(a)}x^{a-1}e^{-\lambda x},\quad x>0,$$

其中 $\Gamma(a)=\int_0^\infty x^{a-1}e^{-x}dx$，$\Gamma(n)=(n-1)!$。课件称$\lambda$为“尺度参数”，但按密度它实际为率参数，讲义统一纠正。

形状为1时就是指数分布。独立Gamma变量在**率相同**时形状相加。形状为正整数$n$时，可表示为$n$个独立同率指数变量之和；非整数形状不能按“非整数个随机变量相加”理解。C3 p22选项C缺少整数与独立同率条件，不能当作一般正确命题。

## 4.7 第$n$次发生时刻与相邻间隔

【课件 C3 p24–32】定义

$$T_n=\text{第}n\text{次事件的发生时刻},\quad T_0=0,\quad X_n=T_n-T_{n-1}.$$

于是 $T_n=X_1+\cdots+X_n$。

课件定理：$X_1,X_2,\ldots$ 相互独立，均服从 $\operatorname{Exp}(\lambda)$。从而

$$T_n\sim\operatorname{Gamma}(n,\lambda),\quad
f_{T_n}(t)=\frac{\lambda^n}{(n-1)!}t^{n-1}e^{-\lambda t},\ t>0.$$

对应 $E(T_n)=n/\lambda$、$D(T_n)=n/\lambda^2$。这两个矩是同一密度的直接推导，不是新增过程模型。

第一个等待时间的推导最直观：$\{T_1>t\}=\{N(t)=0\}$，所以 $P(T_1>t)=e^{-\lambda t}$。一般地，

$$\{T_n>t\}=\{N(t)\le n-1\}.$$

因此等候第$n$次事件超过$t$的概率，也可写成 $e^{-\lambda t}\sum_{k=0}^{n-1}(\lambda t)^k/k!$。

课件给出的第三种描述：若到达间隔独立同分布为指数，所生成的计数过程是Poisson过程。它与“增量为Poisson分布”“小区间条件”相对应。

**汽车例（p27–28）：** 相邻红车间隔率为$\lambda_1$；不分颜色的相邻汽车间隔率为$\lambda_1+\lambda_2+\lambda_3$。这里“相邻”的对象不同，不能混用率。

**仪器振动故障例（p31）：** 第$k$次振动即故障，故障时刻是$T_k$；在$t_0$仍正常的概率为

$$P(T_k>t_0)=e^{-\lambda t_0}\sum_{j=0}^{k-1}\frac{(\lambda t_0)^j}{j!}.$$

## 4.8 课件中的率估计说明

【课件 C3 p26、p30】如果观察$n$个完整独立间隔，累计到第$n$次到达时刻$t_n$，似然为 $\lambda^ne^{-\lambda t_n}$，故极大似然估计 $\hat\lambda=n/t_n$。

若题意是固定观察到时刻$T$、期间共$n$次事件，则应计入最后一次到达后到$T$的无到达区间，似然相应为常数乘 $\lambda^ne^{-\lambda T}$，得到$n/T$。不要把固定观察窗和“第$n$次到达即停止”混为一谈。这里只澄清课件本页设定，不扩展统计推断课程。

由 $T_n$ 的Gamma分布，$2\lambda T_n\sim\chi^2(2n)$。若 $q_\gamma$ 表示该卡方分布的**左尾**$\gamma$分位数，那么基于$n$个完整间隔的置信水平$1-\alpha$区间为

$$\left[\frac{q_{\alpha/2}}{2T_n},\frac{q_{1-\alpha/2}}{2T_n}\right].$$

原课件使用的卡方下标可能是上尾记号，必须先核实分位数约定再代数值，不可只看下标大小。掌握本页结论和使用条件即可。

## 4.9 课件计数例题

【课件 C3 p33–34】8点开始有充分多顾客排队，仅一名服务员连续服务；每次服务时间独立，指数分布均值20分钟。以小时为单位率$\lambda=3$。连续忙碌条件下，到12点离开人数 $N(4)\sim\operatorname{Poi}(12)$，平均12人，恰好9人为 $e^{-12}12^9/9!$。

注意“充分多顾客”保证服务员不空闲；没有该条件，不能无条件把实际离开人数当作此Poisson过程。

## 4.10 给定次数后的发生时刻：均匀与排序

【课件 C3 p35–38】给定 $N(t)=1$，唯一的发生时刻在$(0,t)$均匀分布。对$0\le s\le t$，

$$P(T_1\le s\mid N(t)=1)=\frac{P(N(s)=1)P(N(t)-N(s)=0)}{P(N(t)=1)}=\frac st.$$

更一般地，给定 $N(t)=n$，按时间排好顺序的发生时刻联合密度为

$$f(t_1,\ldots,t_n\mid N(t)=n)=\frac{n!}{t^n},\quad0<t_1<\cdots<t_n<t,$$

其他位置为0。它等于$n$个独立$U(0,t)$变量的**顺序统计量**的联合分布。

**易错点：** 排好顺序的 $T_1,\ldots,T_n$ 不独立，因为必须有先后关系；“独立均匀”描述的是排序之前的抽样点，不是第1、第2次发生时刻各自独立。

## 4.11 题型识别与掌握检查

|题目问法|先选择的对象|
|---|---|
|一段时间内发生几次|Poisson计数或增量|
|下一次还要等多久|指数间隔|
|第$n$次何时发生|Gamma等待时刻|
|给定区间内只发生一次，发生在何处|条件均匀|
|给定共$n$次的有序时刻|均匀顺序统计量|
|多类到达合并、独立筛选|叠加、分解|

需要理解：数字特征的协方差推导、分解的全概率求和、条件均匀的证明。需要会用：增量概率、分解叠加、$T_n$与$X_n$的区别、两类汽车间隔、故障时刻例。视频5说明部分证明因课时限制略去，不能把本讲义补出的所有推导自动变成必须背诵的证明要求。

**同范围检查：** 率为2/小时，求：①半小时恰有1次；②两小时无到达；③第3次到达的平均时刻；④给定两小时内恰有1次，它在前半小时内发生的概率。

教师答案：①$e^{-1}$；②$e^{-4}$；③$3/2$小时；④$1/4$。第④不能直接用未加条件的 $P(N(0.5)\ge1)$。
