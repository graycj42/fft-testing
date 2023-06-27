# Vitis Tutorial

This github is an example vector addition kernel.

## Installation 

```
git clone https://github.com/pyuvaraj37/vitis_tutorial
cd vitis_tutorial
```

## Host 

The host program is located in **src** and takes 1 arguments which is the bitstream to program the connected device. The helper functions for programming the FPGA are in **include/host.h**, which are taken from *Vitis_Accel_Examples*. 

## Kernel 

The Kernel reads from two buffers, which can be either connected to DDR or HBM; depending on the device used. The vectors are then added and written into a third memory buffer. The ports are master AXI interfaces. All the additions are done in parallel using a Vitis HLS pragam **UNROLL**. 

## Synthesizing and Compiling

The Makefile can be used to make both the host and kernel program. The host program is compiled with g++ and OpenCL. The kernel uses Vitis (v++) to compile and depending on software emulation, hardware emulation, and hardware the synthesis time will vary. Before you can run make, the Vitis and XRT versions need to be set.

```
source /tools/Xilinx/Vitis/20XX.X/settings64.sh
source /opt/xilinx/xrt/setup.sh
```

Inside the Makefile will have different configurations such as run, build, or host. When compiling the kernel a TARGET needs to be set, there are three options sw_emu, hw_emu and hw. sw_emu is software emulation and should be used to quickly verify software correctness (Compiles quickly and runs quickly), cannot be used to verify hardware correctness. hw_emu is hardware emulation, and simulates the hardware implmentation to verifty hardware correctness (compiles quickly but runs slow). hw generates a bitstream that programs the device (compiles very slowly, and runs fast). Using a combination of these different TARGETs, we can debug our design. 

The PLATFORM is the specific device that we are developing for. Generally **/opt/xilinx/platforms/** is the location of all the platforms that are installed on the machine. We then have to provide the make file a platform file (.xpfm). 

I would recommend adding the following lines to your **.bashrc**

```
source /tools/Xilinx/Vitis/20XX.X/settings64.sh
source /opt/xilinx/xrt/setup.sh
```

If you will not be using different devices you can add this to your **.bashrc** aswell: 

```
export PLATFORM=/opt/xilinx/platforms/<device>/<device>.xpfm
```


You could also export the TARGET but this will change more often. 

## Make commands

```
make run TARGET=sw_emu/hw_emu/hw PLATFORM=/opt/xilinx/platforms/*device*/*device*.xpfm
```

This command builds both the host and kernel, and then executes the design for the speicifed TARGET. 

```
make build TARGET=sw_emu/hw_emu/hw PLATFORM=/opt/xilinx/platforms/*device*/*device*.xpfm
```
Builds only the kernel program. 

```
make host 
```

Only builds the host program. 

## Config File 
The config file **config.cfg** is 1 of 2 ways to control how the Vitis compiler syntheisizes the kernel to hardware. For example the connectivity of the FPGA design can be specified. In this example, we show the vector memory buffers can be instaniated in HBM or DDR. Other configuration options can be found https://docs.xilinx.com/r/en-US/ug1393-vitis-application-acceleration/v-General-Options. There are many options for profiling, debugging, etc. 

## HLS PRAGMAS 
The other way to control how the kernel is synthesized to hardware is using HLS PRAGMAS (#pragma), these commands are placed within the kernel HLS code. In this example we use only a 2 different types of pragmas. One is an INTERFACE pragma that specify the a, b, and c ports as AXI master ports. The AXI protocol is a standard and can be read about more https://support.xilinx.com/s/article/1053914?language=en_US. The other pragma used is the UNROLL pragma. A common technique in software is to unroll loops to perform more computation per iteration. Thhe UNROLL pragma acheives something similar, by adding more hardware units to perform the computation in parallel. As you might expect this pragma is limited by loop dependencies. The UNROLL pragma can also take in parameters such as factor which specifies the amount of unrolling to perform. For example some extremely large loops cannot be unrolled fully since there will be a physical resource limitation. Try different unroll factors to see the changes in performance. Other HLS pragmas can be found https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/HLS-Pragmas. 

## Tutorial Instructions

We are going to work with the Alveo U50 Data Center Accelerator Card which is installed in Wolverine. The documentation and datasheet can be found here: https://www.xilinx.com/products/boards-and-kits/alveo/u250.html#documentation. 

If you haven't already specified the versions of the tools you are going to be using by executing: 
```
source /tools/Xilinx/Vitis/20XX.X/settings64.sh
source /opt/xilinx/xrt/setup.sh
```
As mentioned above you can add these to .bashrc to automatically run them on log in. 

### Step 1: Software Emulation

Taking a look at the kernel, we are performing a vector addition. The array's a and b are being added into c.

KERNEL
```C++
extern "C" {
    void krnl(data_t* a, data_t* b, data_t* c) {
        #pragma HLS INTERFACE m_axi port = a bundle = gmem0
        #pragma HLS INTERFACE m_axi port = b bundle = gmem1
        #pragma HLS INTERFACE m_axi port = c bundle = gmem0

        for (int i = 0; i < DATA_SIZE; i++) {
            #pragma HLS UNROLL //factor=2
            c[i] = a[i] + b[i];
        }
    }
}
```
HOST
```C++
    for (int i = 0; i < DATA_SIZE; i++) {
        c_sw[i] = a[i] + b[i];
    }
```
We can confirm the addition is providing the correct results in software by running: 
```
make run TARGET=sw_emu PLATFORM=/opt/xilinx/platforms/xilinx_u250_gen3x16_xdma_4_1_202210_1/xilinx_u250_gen3x16_xdma_4_1_202210_1.xpfm 
```
If the the TEST PASSED, we can confirm the values are correct. 

### Step 2: Hardware Emulation
To confirm the hardware works properly we can run a hardware emulation with: 
```
make run TARGET=hw_emu PLATFORM=/opt/xilinx/platforms/xilinx_u250_gen3x16_xdma_4_1_202210_1/xilinx_u250_gen3x16_xdma_4_1_202210_1.xpfm 
```
This may take longer to compiler and run than the sw_emu. A TEST PASSED should be displayed after running. The hardware emulation provides us with an estimate of the hardware resources used. The System estimate can be found in the build folder. 

```
cd _x.hw_emu.xilinx_u250_gen3x16_xdma_4_1_202210_1/reports/link/
vim system_estimate_krnl.link.xtxt
```

The HLS PRAGMA UNROLL can take a few arguments. The one of interest in this tutorial is *factor*. Try uncommenting the **factor=2** in *krnl.cpp* and rerun the hardware emulation. Then try different factors or remove the UNROLL pragma. 

### Step 3: Hardware 
We have only been emulated the hardware until now. 
```
make run TARGET=hw PLATFORM=/opt/xilinx/platforms/xilinx_u250_gen3x16_xdma_4_1_202210_1/xilinx_u250_gen3x16_xdma_4_1_202210_1.xpfm 
```
This command will take a signficantly longer time (Use TMUX if needed). 

Once that completes, look at the system estimate file. Also we can get the final hardware resources used. 

```
cd _x.hw.xilinx_u250_gen3x16_xdma_4_1_202210_1/reports/link/imp/
vim impl_1_full_util_routed.rpt
vim impl_1_hw_bb_locked_timing_summary_routed.rpt
```

These files show the resource utilization and timing reports for the final design which is placed and routed. 

Try different UNROLL factors and try the PIPELINE #pragma to see the differences they make with the resource utilization and performance. 

## More Info & Tips

The build folder will have a report for the synthesized kernel. There you can see the latency and resource utilization of the kernel. 

If you're using git, clean the project before commiting the changes also look for hidden folders **.run** and **.ipcache** these folders will get very large and cause issues pushing your repo. Some bitstreams (.xclbin) files will get too large for Git. 

If you notice a hw compile is taking a long(er) time, try a sw_emu or hw_emu to see if there are any errors. Sometimes when synthezing for hw the software overlooks syntax errors. 

If integrating Verilog/VHDL designs with your HLS design, you can only perform hw_emu and hw runs. 

A powerful tool to use for HW synthesis is TMUX https://www.redhat.com/sysadmin/introduction-tmux-linux. Sometimes HW synthesis can take multiple hours, and maintaining a ssh command for an extended period of time can be difficult. So you can use tmux 
to allow the synthesis continue without an SSH connection. The job can be disowned as well, but TMUX is prefereable as you can reconnect to the tmux session to check the progress easily.  
