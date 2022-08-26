# TempFEM_Parallel

#### 介绍
用TLM算法、DD-TLM算法和Adaptive-TLM计算温度场

#### 软件架构
new_tempFEM：使用单线程直接法矩阵SuperLU，通过DD-TLM的区域分解实现并行化
tempFEM_Paraller，主要是传统TLM和Adaptive_TLM，使用并行直接法矩阵求解器SuperLU_MT

#### 编译环境：
- visual studio2019及以上：https://visualstudio.microsoft.com/zh-hans/
- Windows SDK 10.0以上+MSVC v142生成工具：打开visual studio installer, 勾选Windows 10 SDK和MSVC v142 - VS 2019 C++ x64/x86生成工具
- Qt 5.14.2 + Qt VS Addin：https://download.qt.io/archive/qt/ Qt安装时选择MSVC2019 64bit编译器
