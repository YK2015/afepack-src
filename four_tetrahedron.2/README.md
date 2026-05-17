# four_tetrahedron.2 案例工作区

这个目录承载 `four_tetrahedron.2` 的教学实现与验证材料。

当前已经包含：

- `four_tetrahedron.2.tmp_dof`
- `four_tetrahedron.2.bas_fun`
- `four_tetrahedron.2.bas_fun.c`

这套样例的核心设计是：

- 把 `four_tetrahedron.1` 看成 4 个子四面体上的连续分片一次单元；
- 用现有 7 个一次 nodal basis `lambda_i` 生成二次基函数；
- 顶点基函数使用 `lambda_i(2 lambda_i - 1)`；
- 边中点基函数使用 `4 lambda_i lambda_j`；
- 对宏模板中没有显式编号的 3 条内部微边 `(0,4)`, `(0,5)`, `(0,6)`，
  将其中点 DOF 分别挂到包含它们的宏 face 上。

因此：

- `four_tetrahedron.2.tmp_dof` 里共有 22 条附着记录；
- 总基函数数是 22；
- `2 4 1` 表示在第 4 个宏 face 上挂 `(0,4)` 中点 DOF；
- `2 5 1` 表示在第 5 个宏 face 上挂 `(0,5)` 中点 DOF；
- `2 6 1` 表示在第 6 个宏 face 上挂 `(0,6)` 中点 DOF。

如果需要把该样例编译成共享库，可在当前目录执行：

```bash
cc -shared -fPIC -O2 -o four_tetrahedron.2.bas_fun.so four_tetrahedron.2.bas_fun.c
```

做完编译后，可进一步运行插值验证：

```bash
python3 verify_four_tetrahedron_2.py
```

当前实现已经通过这一步检查，22 个插值点上的基函数矩阵为单位阵。

该案例在手册中对应：

- `template-element-handbook/chapters/10-case-four-tetrahedron-p2.tex`
