# twin_tetrahedron.2 案例工作区

这个目录承载 `twin_tetrahedron.2` 的教学实现与验证材料。

当前已经包含：

- `twin_tetrahedron.2.tmp_dof`
- `twin_tetrahedron.2.bas_fun`
- `twin_tetrahedron.2.bas_fun.c`
- `verify_twin_tetrahedron_2.py`

这套样例的核心设计是：

- 将 `twin_tetrahedron.1` 看成两块子四面体上的连续分片一次单元；
- 用 5 个一次 nodal basis `lambda_i` 生成连续 piecewise P2 基函数；
- 7 条显式宏边上的中点 DOF 继续挂在 edge entity 上；
- 额外 2 条内部微边 `(0,2)` 和 `(2,4)` 的中点 DOF 分别挂到包含它们的宏 face 上。

因此：

- `twin_tetrahedron.2.tmp_dof` 里共有 14 条附着记录；
- 总基函数数是 14；
- `2 0 1` 表示在第 0 个宏 face 上挂 `(2,4)` 中点 DOF；
- `2 3 1` 表示在第 3 个宏 face 上挂 `(0,2)` 中点 DOF。

如需编译共享库，可在当前目录执行：

```bash
cc -shared -fPIC -O2 -o twin_tetrahedron.2.bas_fun.so twin_tetrahedron.2.bas_fun.c
```

编译后可做插值验证：

```bash
python3 verify_twin_tetrahedron_2.py
```

该案例在手册中对应：

- `template-element-handbook/chapters/11-case-twin-tetrahedron-p2.tex`
