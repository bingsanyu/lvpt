# evaluation.R 与 evaluation.py

## *待解决ToDo*

When run

```R
models$model <- map(models$model, add_cell_waypoints)
```

get the error

```R
Error in `map()`:
ℹ In index: 1.
Caused by error:
! is_wrapper_with_trajectory(trajectory = trajectory) is not TRUE
```

## 新建Python虚拟环境

```Bash
conda create -n lvpt python=3.8 scanpy scvelo
```

## 修改代码

在`evaluation.R`开头添加以下内容

```R
# 开启代理
if (!require(r.proxy)) {
  install.packages("r.proxy")
}
r.proxy::proxy()
# 安装R包
if (!require(devtools)) {
  install.packages("devtools")
}
if (!require(dyno)) {
  devtools::install_github("dynverse/dyno")
}
if (!require(dyneval)) {
  devtools::install_github("dynverse/dyneval")
}
# 创建用于存放输出数据的文件夹
simulated_result_dir = "simulated_result"
simulated_data_dir = "simulated_data"
if (!dir.exists(simulated_result_dir)){
  dir.create(simulated_result_dir)
}
if (!dir.exists(simulated_data_dir)){
  dir.create(simulated_data_dir)
}
# 指定Python环境（根据你的anaconda安装位置进行修改）
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "D:/anaconda/envs/lvpt/python.exe")
use_python("D:/anaconda/envs/lvpt/python.exe")
py_config()
```

在`res = c()`后面一行添加

```R
source_python("./evaluation.py")
```

## 运行代码

在Rstudio中运行`evaluation.R`。
