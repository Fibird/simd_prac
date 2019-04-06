# SIMD 练习

这里的几个代码都是SIMD学习的几个练习，包括数据库的一些基本操作的SIMD实现。

# TPCH Q1的SIMD实现

代码详见vectorized_process.c和vectorized_process_v2.c。

# Nest Loop Join的优化

## 基于SIMD的优化

代码详见simd_nlj.c、simd_nlj_v2.c。

## 基于划分（partition）的优化

代码详见part_nlj.c、simd_part_nlj.c。

# 编译方法

在编译指令后使用`-mavx`或者`-mavx2`来对这些程序进行编译

```
gcc -mavx -O1 <source>.c -o <exe_name>
```

# 学习资料

- [Intel Intrinsics Guide](https://software.intel.com/sites/landingpage/IntrinsicsGuide/).
