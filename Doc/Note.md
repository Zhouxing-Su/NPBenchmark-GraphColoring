先贪心对未覆盖的节点找最大独立集 减少无用的初始集合数

**图着色局部搜索 桶队列 解禁忌 固定最大团**
  **最大团选择在节点数相同的情况下度最大的？**

**不能覆盖所有节点时如何生成新独立集？**
  生成独立集时考虑节点的度？

**优度优先搜索**
  两个节点在同一集合/不在同一集合
  优度=潜在剩余冲突数（下界估计或向前看树搜索）

subgraph的键改为子图节点集或者应尽量包含的节点集（与IRP里面类似）？

图着色已知任意两点间的颜色编号大小关系的情况下求解是否是P问题？
  对编号大小关系分支限界？
    **最大团内的节点**
图着色已知任意两点间的颜色是否相同的情况下求解是否是P问题？
  仍然是一个图着色问题，但约束更紧

**禁忌边不能重新产生冲突**

解禁忌编码为两个节点是否在同一颜色的布尔向量
  等价于记录冲突边编号
  太长？
    只记录 i < j 的项（下三角矩阵）
    记录非零元的下标数组

每个节点找个最大独立集 去掉重复的 算集合覆盖 对没覆盖的节点分支计算次优独立集
  枚举次优独立集算法
    分支排除最大独立集中的每个节点 优度优先搜索
  针对未覆盖的原因构造独立集（不一定是对未覆盖的点算次优独立集）
    提升未覆盖节点的权重
      计算独立集时
      未覆盖的惩罚
    包含一个未覆盖点和一个已覆盖点的最大独立集
    所有独立集中缺乏的节点组合
      某个节点只在很少的独立集中

图着色交叉算符精确计算最大平衡继承/P-Center交叉算符优化原始目标平衡继承
进化计算与列生成的平衡点，甚至可以不需要局部搜索
甚至不用保存完整解，直接变成一大堆独立集组合，用贪心算法，不一定挑2K个集合（每个解各K个），可以挑更多
  完整解里面每个颜色的节点集合其实是可能有冲突的，不一定是独立集
  江华精确算法/向前看树搜索启发式搜索最大独立集
直接优化冲突数或者分两阶段先尽量多覆盖节点再优化冲突数（原本集合内的冲突+未覆盖节点的冲突）

精确算法
  向前看树搜索？DSATUR
  每个节点相邻未染色节点最大团快速依次染色
  从松弛初始解出发尝试变换当前解
    去重？消除对称性？
  蒙特卡洛树搜索？
  分支策略: 优先选择难满足的分支
    剩余可用颜色数最少
    未染色相邻节点数最少
    预计不可用颜色总数 = (已用颜色数 / 已染色相邻节点数) * 总相邻节点数
  **将节点按颜色展开, 相邻节点不同颜色之间有虚拟边, 不相邻节点任意颜色之间有边, 求最大团**
    稀疏图上求独立集/稠密图上求最大团
    消除对称性的冗余约束消除无用虚拟节点和虚拟边
      原始图上的最大团直接固定
      颜色编号与节点编号的关系

约束的 price 和 slack 有什么实际意义？
图着色判定版本中多出的独立集数量约束直接忽略其 price 求子问题是否会有问题？
