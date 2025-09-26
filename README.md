# SAC Adaptive Modulation Physical Layer (MATLAB)

## 项目简介 / Project Description
本仓库提供基于 MATLAB 的 OTFS 物理层仿真，用于 SAC（Soft Actor-Critic）自适应调制系统。该模块作为强化学习环境（Environment），支持延迟-多普勒信道建模、OTFS 调制/解调、消息传递检测以及 BER 计算。

This repository provides MATLAB code for simulating the OTFS physical layer for SAC (Soft Actor-Critic) adaptive modulation. This module serves as a reinforcement learning environment, supporting delay-Doppler channel modeling, OTFS modulation/demodulation, message passing detection, and BER calculation.

---

## 功能 / Features
- OTFS 调制与解调 / OTFS modulation & demodulation
- 延迟-多普勒信道建模 / Delay-Doppler channel modeling
- 消息传递检测算法 / Message passing detection algorithm
- 输出比特误码率（BER）作为 RL 奖励 / Output Bit Error Rate (BER) as RL reward
- 可根据 SNR 变化生成环境状态 / Generate environment states based on SNR
- 支持不同调制阶数动态映射 / Support dynamic modulation order mapping

---

## 文件结构 / File Structure
SAC_OTFS_PHY/
│-- main_otfs_simulation.m % 主物理层仿真代码
│-- otfs_data.mat % 信道数据文件
│-- README.md % 项目说明文档
