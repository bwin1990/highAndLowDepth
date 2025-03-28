#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DNA寡核苷酸自我互补结构分析工具 - 图形用户界面
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
from dna_analyzer import DNAComplementAnalyzer, ASusceptibilityAnalyzer
from typing import List, Dict, Any
import json
import platform
import matplotlib.font_manager as fm
import math

# 配置matplotlib支持中文
if platform.system() == 'Windows':
    # Windows系统字体配置
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'SimSun', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    # 尝试加载微软雅黑字体
    font_path = 'C:/Windows/Fonts/msyh.ttc'  # 微软雅黑字体路径
    if os.path.exists(font_path):
        font_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = font_prop.get_name()
else:
    # 其他系统字体配置
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'WenQuanYi Micro Hei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False


def read_sequences_from_file(file_path: str) -> List[str]:
    """从文件中读取DNA序列"""
    sequences = []
    
    # 检查文件扩展名
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.fasta' or ext == '.fa':
        # 读取FASTA格式
        with open(file_path, 'r') as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = ""
                else:
                    current_seq += line
            if current_seq:
                sequences.append(current_seq)
    
    elif ext == '.txt':
        # 每行一个序列
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    sequences.append(line)
    
    elif ext == '.csv':
        # 试图读取CSV文件，假设第一列是序列
        df = pd.read_csv(file_path)
        if len(df.columns) > 0:
            sequences = df.iloc[:, 0].tolist()
    
    elif ext == '.tsv':
        # 试图读取TSV文件，假设第一列是序列
        df = pd.read_csv(file_path, sep='\t')
        if len(df.columns) > 0:
            sequences = df.iloc[:, 0].tolist()
    
    else:
        raise ValueError(f"不支持的文件格式: {ext}")
    
    return sequences

class OligoAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA寡核苷酸自我互补结构分析工具")
        self.root.geometry("900x700")
        self.root.minsize(800, 600)
        
        # 添加窗口关闭事件处理
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # 创建分析器实例
        self.analyzer = DNAComplementAnalyzer()
        self.a_analyzer = ASusceptibilityAnalyzer()
        
        # 创建主框架
        self.main_frame = ttk.Notebook(root)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 创建单序列分析标签页
        self.single_frame = ttk.Frame(self.main_frame)
        self.main_frame.add(self.single_frame, text="单序列分析")
        self._setup_single_analysis_tab()
        
        # 创建组比较标签页
        self.compare_frame = ttk.Frame(self.main_frame)
        self.main_frame.add(self.compare_frame, text="组比较分析")
        self._setup_compare_analysis_tab()
        
        # 创建A碱基susceptibility分析标签页
        self.a_suscept_frame = ttk.Frame(self.main_frame)
        self.main_frame.add(self.a_suscept_frame, text="A碱基Susceptibility分析")
        self._setup_a_suscept_tab()
        
        # 创建状态栏
        self.status_var = tk.StringVar()
        self.status_var.set("就绪")
        self.status_bar = ttk.Label(root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
    def _on_closing(self):
        """处理窗口关闭事件"""
        # 清理资源
        plt.close('all')  # 关闭所有matplotlib图形
        
        # 销毁根窗口并退出程序
        self.root.destroy()
        self.root.quit()  # 确保退出主循环
        
    def _setup_single_analysis_tab(self):
        """设置单序列分析标签页"""
        # 创建左侧输入区域
        input_frame = ttk.LabelFrame(self.single_frame, text="输入")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 序列输入
        ttk.Label(input_frame, text="DNA序列:").pack(anchor=tk.W, padx=5, pady=2)
        self.sequence_text = scrolledtext.ScrolledText(input_frame, height=10)
        self.sequence_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 参数设置
        params_frame = ttk.Frame(input_frame)
        params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(params_frame, text="Flank长度:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.flank_var = tk.IntVar(value=20)
        ttk.Spinbox(params_frame, from_=5, to=50, textvariable=self.flank_var, width=10).grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="窗口大小:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.window_var = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=5, to=30, textvariable=self.window_var, width=10).grid(row=1, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="最小匹配:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        self.min_match_var = tk.IntVar(value=4)
        ttk.Spinbox(params_frame, from_=3, to=15, textvariable=self.min_match_var, width=10).grid(row=2, column=1, padx=5, pady=2)
        
        # 分析模式选择
        ttk.Label(params_frame, text="分析模式:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.analysis_mode_var = tk.StringVar(value="flank")
        analysis_mode_combo = ttk.Combobox(params_frame, textvariable=self.analysis_mode_var, width=10)
        analysis_mode_combo['values'] = ('flank', 'global')
        analysis_mode_combo['state'] = 'readonly'
        analysis_mode_combo.grid(row=3, column=1, padx=5, pady=2)
        analysis_mode_combo.bind('<<ComboboxSelected>>', self._on_analysis_mode_change)
        
        # Flank分析模式选择
        ttk.Label(params_frame, text="Flank分析模式:").grid(row=4, column=0, sticky=tk.W, padx=5, pady=2)
        self.flank_mode_var = tk.StringVar(value="both")
        self.flank_mode_combo = ttk.Combobox(params_frame, textvariable=self.flank_mode_var, width=10)
        self.flank_mode_combo['values'] = ('both', '5p', '3p')
        self.flank_mode_combo['state'] = 'readonly'
        self.flank_mode_combo.grid(row=4, column=1, padx=5, pady=2)
        
        # GC-only选项（只在全局模式下显示）
        self.gc_only_var = tk.BooleanVar(value=False)
        self.gc_only_check = ttk.Checkbutton(params_frame, text="仅分析GC序列", variable=self.gc_only_var)
        self.gc_only_check.grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=5, pady=2)
        self.gc_only_check.state(['disabled'])  # 初始时禁用
        
        # 添加模式说明
        mode_frame = ttk.Frame(input_frame)
        mode_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(mode_frame, text="分析模式: flank=分析flank与内部序列互补, global=全局序列互补分析").pack(anchor=tk.W)
        ttk.Label(mode_frame, text="Flank模式: both=两端flank, 5p=仅5'端flank, 3p=仅3'端flank").pack(anchor=tk.W)
        
        # 按钮区域
        button_frame = ttk.Frame(input_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=10)
        
        ttk.Button(button_frame, text="分析序列", command=self._analyze_single_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除", command=lambda: self.sequence_text.delete(1.0, tk.END)).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存结果", command=self._save_single_results).pack(side=tk.LEFT, padx=5)
        
        # 创建右侧结果区域
        result_frame = ttk.LabelFrame(self.single_frame, text="结果")
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 结果文本区域
        self.result_text = scrolledtext.ScrolledText(result_frame, height=10)
        self.result_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 图表区域
        self.figure_frame = ttk.Frame(result_frame)
        self.figure_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 初始化结果存储
        self.single_results = None
        self.single_score = None
        
    def _on_analysis_mode_change(self, event):
        """分析模式改变时的处理函数"""
        if self.analysis_mode_var.get() == 'flank':
            # 启用Flank模式选项，禁用GC-only选项
            self.flank_mode_combo.state(['!disabled'])
            self.gc_only_check.state(['disabled'])
            self.gc_only_var.set(False)
        else:  # global模式
            # 禁用Flank模式选项，启用GC-only选项
            self.flank_mode_combo.state(['disabled'])
            self.gc_only_check.state(['!disabled'])
        
    def _setup_compare_analysis_tab(self):
        """设置组比较分析标签页"""
        # 创建左侧输入区域
        input_frame = ttk.LabelFrame(self.compare_frame, text="输入")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 高效率组文件选择
        file_frame1 = ttk.Frame(input_frame)
        file_frame1.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(file_frame1, text="高效率组文件:").pack(side=tk.LEFT, padx=5)
        self.high_eff_path_var = tk.StringVar()
        ttk.Entry(file_frame1, textvariable=self.high_eff_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(file_frame1, text="浏览...", command=lambda: self._browse_file(self.high_eff_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 低效率组文件选择
        file_frame2 = ttk.Frame(input_frame)
        file_frame2.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(file_frame2, text="低效率组文件:").pack(side=tk.LEFT, padx=5)
        self.low_eff_path_var = tk.StringVar()
        ttk.Entry(file_frame2, textvariable=self.low_eff_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(file_frame2, text="浏览...", command=lambda: self._browse_file(self.low_eff_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 参数设置
        params_frame = ttk.Frame(input_frame)
        params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(params_frame, text="Flank长度:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_flank_var = tk.IntVar(value=20)
        ttk.Spinbox(params_frame, from_=5, to=50, textvariable=self.comp_flank_var, width=10).grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="窗口大小:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_window_var = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=5, to=30, textvariable=self.comp_window_var, width=10).grid(row=1, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="最小匹配:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_min_match_var = tk.IntVar(value=4)
        ttk.Spinbox(params_frame, from_=3, to=15, textvariable=self.comp_min_match_var, width=10).grid(row=2, column=1, padx=5, pady=2)
        
        # 分析模式选择
        ttk.Label(params_frame, text="分析模式:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_analysis_mode_var = tk.StringVar(value="flank")
        comp_analysis_mode_combo = ttk.Combobox(params_frame, textvariable=self.comp_analysis_mode_var, width=10)
        comp_analysis_mode_combo['values'] = ('flank', 'global')
        comp_analysis_mode_combo['state'] = 'readonly'
        comp_analysis_mode_combo.grid(row=3, column=1, padx=5, pady=2)
        comp_analysis_mode_combo.bind('<<ComboboxSelected>>', self._on_comp_analysis_mode_change)
        
        # Flank分析模式选择
        ttk.Label(params_frame, text="Flank分析模式:").grid(row=4, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_flank_mode_var = tk.StringVar(value="both")
        self.comp_flank_mode_combo = ttk.Combobox(params_frame, textvariable=self.comp_flank_mode_var, width=10)
        self.comp_flank_mode_combo['values'] = ('both', '5p', '3p')
        self.comp_flank_mode_combo['state'] = 'readonly'
        self.comp_flank_mode_combo.grid(row=4, column=1, padx=5, pady=2)
        
        # GC-only选项（只在全局模式下显示）
        self.comp_gc_only_var = tk.BooleanVar(value=False)
        self.comp_gc_only_check = ttk.Checkbutton(params_frame, text="仅分析GC序列", variable=self.comp_gc_only_var)
        self.comp_gc_only_check.grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=5, pady=2)
        self.comp_gc_only_check.state(['disabled'])  # 初始时禁用
        
        # 图表类型选择
        ttk.Label(params_frame, text="图表类型:").grid(row=6, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_plot_type_var = tk.StringVar(value="boxviolin")
        comp_plot_type_combo = ttk.Combobox(params_frame, textvariable=self.comp_plot_type_var, width=10)
        comp_plot_type_combo['values'] = ('boxviolin', 'histogram', 'scatter', 'boxplot', 'violinplot', 'swarmplot', 'stripplot')
        comp_plot_type_combo['state'] = 'readonly'
        comp_plot_type_combo.grid(row=6, column=1, padx=5, pady=2)
        
        # 添加图表类型说明
        chart_frame = ttk.Frame(input_frame)
        chart_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(chart_frame, text="分析模式: flank=分析flank与内部序列互补, global=全局序列互补分析").pack(anchor=tk.W)
        ttk.Label(chart_frame, text="图表类型: boxviolin=箱线图+小提琴图, histogram=直方图, scatter=散点图").pack(anchor=tk.W)
        
        # 添加模式说明
        mode_frame = ttk.Frame(input_frame)
        mode_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(mode_frame, text="Flank模式: both=两端flank, 5p=仅5'端flank, 3p=仅3'端flank").pack(anchor=tk.W)
        
        # 输出设置
        output_frame = ttk.Frame(input_frame)
        output_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(output_frame, text="结果输出文件:").pack(side=tk.LEFT, padx=5)
        self.output_path_var = tk.StringVar()
        ttk.Entry(output_frame, textvariable=self.output_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(output_frame, text="浏览...", command=lambda: self._save_file_dialog(self.output_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 按钮区域
        button_frame = ttk.Frame(input_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=10)
        
        ttk.Button(button_frame, text="比较分析", command=self._compare_groups).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除", command=self._clear_compare_inputs).pack(side=tk.LEFT, padx=5)
        
        # 创建右侧结果区域
        result_frame = ttk.LabelFrame(self.compare_frame, text="结果")
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 结果文本区域
        self.comp_result_text = scrolledtext.ScrolledText(result_frame, height=10)
        self.comp_result_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 图表区域
        self.comp_figure_frame = ttk.Frame(result_frame)
        self.comp_figure_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 初始化结果存储
        self.comparison_results = None
    
    def _on_comp_analysis_mode_change(self, event):
        """组比较分析模式改变时的处理函数"""
        if self.comp_analysis_mode_var.get() == 'flank':
            # 启用Flank模式选项，禁用GC-only选项
            self.comp_flank_mode_combo.state(['!disabled'])
            self.comp_gc_only_check.state(['disabled'])
            self.comp_gc_only_var.set(False)
        else:  # global模式
            # 禁用Flank模式选项，启用GC-only选项
            self.comp_flank_mode_combo.state(['disabled'])
            self.comp_gc_only_check.state(['!disabled'])
        
    def _setup_a_suscept_tab(self):
        """设置A碱基Susceptibility分析标签页"""
        # 创建嵌套的标签页控件
        a_suscept_notebook = ttk.Notebook(self.a_suscept_frame)
        a_suscept_notebook.pack(fill=tk.BOTH, expand=True)
        
        # 创建单序列分析标签页
        self.a_single_frame = ttk.Frame(a_suscept_notebook)
        a_suscept_notebook.add(self.a_single_frame, text="单序列分析")
        
        # 创建批量比较标签页
        self.a_compare_frame = ttk.Frame(a_suscept_notebook)
        a_suscept_notebook.add(self.a_compare_frame, text="批量比较")
        
        # 设置单序列分析区域
        self._setup_a_single_analysis()
        
        # 设置批量比较区域
        self._setup_a_batch_comparison()
    
    def _setup_a_single_analysis(self):
        """设置A碱基单序列分析区域"""
        # 创建左侧输入区域
        input_frame = ttk.LabelFrame(self.a_single_frame, text="输入")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 序列输入
        ttk.Label(input_frame, text="DNA序列:").pack(anchor=tk.W, padx=5, pady=2)
        self.a_sequence_text = scrolledtext.ScrolledText(input_frame, height=10)
        self.a_sequence_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 按钮区域
        button_frame = ttk.Frame(input_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=10)
        
        ttk.Button(button_frame, text="分析序列", command=self._analyze_a_single_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除", command=lambda: self.a_sequence_text.delete(1.0, tk.END)).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="保存结果", command=self._save_a_single_results).pack(side=tk.LEFT, padx=5)
        
        # 创建右侧结果区域
        result_frame = ttk.LabelFrame(self.a_single_frame, text="结果")
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 结果文本区域
        self.a_result_text = scrolledtext.ScrolledText(result_frame, height=10)
        self.a_result_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 图表区域
        self.a_figure_frame = ttk.Frame(result_frame)
        self.a_figure_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 初始化结果存储
        self.a_single_score = None
        self.a_single_sequence = None
        
    def _setup_a_batch_comparison(self):
        """设置A碱基批量比较区域"""
        # 创建左侧输入区域
        input_frame = ttk.LabelFrame(self.a_compare_frame, text="输入")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 高效率组文件选择
        file_frame1 = ttk.Frame(input_frame)
        file_frame1.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(file_frame1, text="高效率组文件:").pack(side=tk.LEFT, padx=5)
        self.a_high_eff_path_var = tk.StringVar()
        ttk.Entry(file_frame1, textvariable=self.a_high_eff_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(file_frame1, text="浏览...", command=lambda: self._browse_file(self.a_high_eff_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 低效率组文件选择
        file_frame2 = ttk.Frame(input_frame)
        file_frame2.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(file_frame2, text="低效率组文件:").pack(side=tk.LEFT, padx=5)
        self.a_low_eff_path_var = tk.StringVar()
        ttk.Entry(file_frame2, textvariable=self.a_low_eff_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(file_frame2, text="浏览...", command=lambda: self._browse_file(self.a_low_eff_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 图表类型选择
        params_frame = ttk.Frame(input_frame)
        params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(params_frame, text="图表类型:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.a_plot_type_var = tk.StringVar(value="boxviolin")
        plot_type_combo = ttk.Combobox(params_frame, textvariable=self.a_plot_type_var, width=10)
        plot_type_combo['values'] = ('boxviolin', 'histogram', 'scatter', 'boxplot', 'violinplot', 'swarmplot', 'stripplot')
        plot_type_combo['state'] = 'readonly'
        plot_type_combo.grid(row=0, column=1, padx=5, pady=2)
        
        # 添加图表类型说明
        chart_frame = ttk.Frame(input_frame)
        chart_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(chart_frame, text="图表类型: boxviolin=箱线图+小提琴图, histogram=直方图, scatter=散点图").pack(anchor=tk.W)
        
        # 输出设置
        output_frame = ttk.Frame(input_frame)
        output_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(output_frame, text="结果输出文件:").pack(side=tk.LEFT, padx=5)
        self.a_output_path_var = tk.StringVar()
        ttk.Entry(output_frame, textvariable=self.a_output_path_var, width=30).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(output_frame, text="浏览...", command=lambda: self._save_file_dialog(self.a_output_path_var)).pack(side=tk.LEFT, padx=5)
        
        # 按钮区域
        button_frame = ttk.Frame(input_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=10)
        
        ttk.Button(button_frame, text="比较分析", command=self._compare_a_susceptibility).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除", command=self._clear_a_compare_inputs).pack(side=tk.LEFT, padx=5)
        
        # 创建右侧结果区域
        result_frame = ttk.LabelFrame(self.a_compare_frame, text="结果")
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 结果文本区域
        self.a_comp_result_text = scrolledtext.ScrolledText(result_frame, height=10)
        self.a_comp_result_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 图表区域
        self.a_comp_figure_frame = ttk.Frame(result_frame)
        self.a_comp_figure_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # 初始化结果存储
        self.a_comparison_results = None
    
    def _browse_file(self, path_var):
        """打开文件浏览对话框"""
        filetypes = [
            ("所有支持的文件", "*.fasta;*.fa;*.txt;*.csv;*.tsv"),
            ("FASTA文件", "*.fasta;*.fa"),
            ("文本文件", "*.txt"),
            ("CSV文件", "*.csv"),
            ("TSV文件", "*.tsv"),
            ("所有文件", "*.*")
        ]
        filename = filedialog.askopenfilename(filetypes=filetypes)
        if filename:
            path_var.set(filename)
            
    def _save_file_dialog(self, path_var):
        """打开保存文件对话框"""
        filetypes = [
            ("Excel文件", "*.xlsx"),
            ("CSV文件", "*.csv"),
            ("JSON文件", "*.json"),
            ("文本文件", "*.txt"),
            ("所有文件", "*.*")
        ]
        filename = filedialog.asksaveasfilename(filetypes=filetypes)
        if filename:
            path_var.set(filename)
            
    def _analyze_single_sequence(self):
        """分析单个序列"""
        sequence = self.sequence_text.get(1.0, tk.END).strip()
        if not sequence:
            messagebox.showerror("错误", "请输入DNA序列")
            return
            
        try:
            self.status_var.set("正在分析序列...")
            self.root.update_idletasks()
            
            # 更新分析器参数
            self.analyzer.flank_length = self.flank_var.get()
            self.analyzer.window_size = self.window_var.get()
            self.analyzer.min_match = self.min_match_var.get()
            self.analyzer.flank_mode = self.flank_mode_var.get()
            
            # 获取分析模式
            analysis_mode = self.analysis_mode_var.get()
            
            # 根据模式分析序列
            if analysis_mode == 'flank':
                self.single_results, self.single_score = self.analyzer.analyze_sequence(sequence)
            else:  # global模式
                gc_only = self.gc_only_var.get()
                self.single_results, self.single_score = self.analyzer.analyze_sequence_global(sequence, gc_only)
            
            # 显示结果
            self._display_single_results(sequence, analysis_mode)
            
            self.status_var.set("分析完成")
        except Exception as e:
            messagebox.showerror("错误", str(e))
            self.status_var.set("分析失败")
            
    def _display_single_results(self, sequence, mode='flank'):
        """显示单序列分析结果"""
        # 清除结果区域
        self.result_text.delete(1.0, tk.END)
        
        # 显示总评分
        self.result_text.insert(tk.END, f"总互补评分: {self.single_score:.2f}\n")
        self.result_text.insert(tk.END, f"发现 {len(self.single_results)} 个潜在互补区域\n")
        
        if mode == 'flank':
            self.result_text.insert(tk.END, f"分析模式: flank (Flank与内部序列互补)\n")
            self.result_text.insert(tk.END, f"Flank分析模式: {self.analyzer.flank_mode}\n\n")
            
            # 显示flank分析详细结果
            for i, result in enumerate(self.single_results):
                self.result_text.insert(tk.END, f"互补区域 {i+1}:\n")
                self.result_text.insert(tk.END, f"Flank: {result['flank']}, 位置: {result['flank_pos']}\n")
                self.result_text.insert(tk.END, f"内部位置: {result['internal_pos']}\n")
                self.result_text.insert(tk.END, f"Flank序列: {result['flank_seq']}\n")
                self.result_text.insert(tk.END, f"内部序列: {result['internal_seq']}\n")
                self.result_text.insert(tk.END, f"互补序列: {result['complementary_seq']}\n")
                self.result_text.insert(tk.END, f"匹配长度: {result['match_length']} bp\n")
                self.result_text.insert(tk.END, f"评分: {result['score']:.2f}\n")
                self.result_text.insert(tk.END, f"预估自由能: {result['free_energy']:.2f} kcal/mol\n\n")
        else:  # global模式
            gc_only = self.gc_only_var.get()
            self.result_text.insert(tk.END, f"分析模式: global (全局序列互补分析)\n")
            if gc_only:
                self.result_text.insert(tk.END, f"仅分析GC序列: 是\n\n")
            else:
                self.result_text.insert(tk.END, f"仅分析GC序列: 否\n\n")
                
            # 显示全局分析详细结果
            for i, result in enumerate(self.single_results):
                self.result_text.insert(tk.END, f"互补区域 {i+1}:\n")
                self.result_text.insert(tk.END, f"区域1位置: {result['region1_pos']}\n")
                self.result_text.insert(tk.END, f"区域2位置: {result['region2_pos']}\n")
                self.result_text.insert(tk.END, f"区域1序列: {result['region1_seq']}\n")
                self.result_text.insert(tk.END, f"区域2序列: {result['region2_seq']}\n")
                self.result_text.insert(tk.END, f"互补序列: {result['complementary_seq']}\n")
                self.result_text.insert(tk.END, f"匹配长度: {result['match_length']} bp\n")
                self.result_text.insert(tk.END, f"评分: {result['score']:.2f}\n")
                self.result_text.insert(tk.END, f"预估自由能: {result['free_energy']:.2f} kcal/mol\n\n")
            
        # 显示可视化
        for widget in self.figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        fig = plt.figure(figsize=(6, 4))
        
        if mode == 'flank':
            self.analyzer.visualize_structure(sequence, self.single_results)
        else:  # global模式
            gc_only = self.gc_only_var.get()
            self.analyzer.visualize_global_structure(sequence, self.single_results, gc_only)
        
        canvas = FigureCanvasTkAgg(fig, master=self.figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def _save_single_results(self):
        """保存单序列分析结果"""
        if not self.single_results:
            messagebox.showerror("错误", "没有可保存的结果")
            return
            
        filetypes = [
            ("文本文件", "*.txt"),
            ("JSON文件", "*.json"),
            ("所有文件", "*.*")
        ]
        filename = filedialog.asksaveasfilename(filetypes=filetypes, defaultextension=".txt")
        if not filename:
            return
            
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(f"DNA自我互补结构分析结果\n")
                f.write("=" * 50 + "\n\n")
                
                # 获取分析模式
                mode = self.analysis_mode_var.get()
                
                f.write(f"总互补评分: {self.single_score:.2f}\n")
                f.write(f"发现 {len(self.single_results)} 个潜在互补区域\n")
                
                if mode == 'flank':
                    f.write(f"分析模式: flank (Flank与内部序列互补)\n")
                    f.write(f"Flank分析模式: {self.analyzer.flank_mode}\n\n")
                    
                    for i, result in enumerate(self.single_results):
                        f.write(f"互补区域 {i+1}:\n")
                        f.write(f"Flank: {result['flank']}, 位置: {result['flank_pos']}\n")
                        f.write(f"内部位置: {result['internal_pos']}\n")
                        f.write(f"Flank序列: {result['flank_seq']}\n")
                        f.write(f"内部序列: {result['internal_seq']}\n")
                        f.write(f"互补序列: {result['complementary_seq']}\n")
                        f.write(f"匹配长度: {result['match_length']} bp\n")
                        f.write(f"评分: {result['score']:.2f}\n")
                        f.write(f"预估自由能: {result['free_energy']:.2f} kcal/mol\n\n")
                else:  # global模式
                    gc_only = self.gc_only_var.get()
                    f.write(f"分析模式: global (全局序列互补分析)\n")
                    f.write(f"仅分析GC序列: {'是' if gc_only else '否'}\n\n")
                    
                    for i, result in enumerate(self.single_results):
                        f.write(f"互补区域 {i+1}:\n")
                        f.write(f"区域1位置: {result['region1_pos']}\n")
                        f.write(f"区域2位置: {result['region2_pos']}\n")
                        f.write(f"区域1序列: {result['region1_seq']}\n")
                        f.write(f"区域2序列: {result['region2_seq']}\n")
                        f.write(f"互补序列: {result['complementary_seq']}\n")
                        f.write(f"匹配长度: {result['match_length']} bp\n")
                        f.write(f"评分: {result['score']:.2f}\n")
                        f.write(f"预估自由能: {result['free_energy']:.2f} kcal/mol\n\n")
                    
            # 保存图像
            img_filename = os.path.splitext(filename)[0] + ".png"
            fig = plt.figure(figsize=(8, 6))
            
            sequence = self.sequence_text.get(1.0, tk.END).strip()
            if mode == 'flank':
                self.analyzer.visualize_structure(sequence, self.single_results, save_path=img_filename)
            else:  # global模式
                gc_only = self.gc_only_var.get()
                self.analyzer.visualize_global_structure(sequence, self.single_results, gc_only=gc_only, save_path=img_filename)
                
            plt.close(fig)
            
            self.status_var.set(f"结果已保存至 {filename}")
            messagebox.showinfo("保存成功", f"结果已保存至:\n{filename}\n图像已保存至:\n{img_filename}")
        except Exception as e:
            messagebox.showerror("保存失败", str(e))
            
    def _compare_groups(self):
        """比较两组序列"""
        high_eff_file = self.high_eff_path_var.get()
        low_eff_file = self.low_eff_path_var.get()
        
        if not high_eff_file or not low_eff_file:
            messagebox.showerror("错误", "请选择高效率和低效率序列文件")
            return
            
        if not os.path.exists(high_eff_file) or not os.path.exists(low_eff_file):
            messagebox.showerror("错误", "文件不存在")
            return
            
        try:
            self.status_var.set("正在读取序列...")
            self.root.update_idletasks()
            
            # 读取序列
            high_eff_seqs = read_sequences_from_file(high_eff_file)
            low_eff_seqs = read_sequences_from_file(low_eff_file)
            
            if not high_eff_seqs or not low_eff_seqs:
                messagebox.showerror("错误", "无法从文件中读取序列")
                self.status_var.set("读取序列失败")
                return
                
            self.status_var.set(f"正在分析 {len(high_eff_seqs)} 条高效率序列和 {len(low_eff_seqs)} 条低效率序列...")
            self.root.update_idletasks()
            
            # 更新分析器参数
            self.analyzer.flank_length = self.comp_flank_var.get()
            self.analyzer.window_size = self.comp_window_var.get()
            self.analyzer.min_match = self.comp_min_match_var.get()
            self.analyzer.flank_mode = self.comp_flank_mode_var.get()
            
            # 获取分析模式
            mode = self.comp_analysis_mode_var.get()
            gc_only = False
            
            # 如果是全局模式，获取GC-only参数
            if mode == 'global':
                gc_only = self.comp_gc_only_var.get()
            
            # 比较两组
            self.comparison_results = self.analyzer.compare_groups(
                high_eff_seqs, 
                low_eff_seqs,
                mode=mode,
                gc_only=gc_only
            )
            
            # 显示结果
            self._display_comparison_results()
            
            # 保存结果
            output_file = self.output_path_var.get()
            if output_file:
                save_results_to_file(self.comparison_results, output_file)
                self.status_var.set(f"结果已保存至 {output_file}")
            else:
                self.status_var.set("分析完成")
                
        except Exception as e:
            messagebox.showerror("错误", str(e))
            self.status_var.set("分析失败")
            
    def _display_comparison_results(self):
        """显示组比较分析结果"""
        # 清除结果区域
        self.comp_result_text.delete(1.0, tk.END)
        
        # 显示分析模式信息
        mode = self.comp_analysis_mode_var.get()
        gc_only = self.comp_gc_only_var.get() if mode == 'global' else False
        
        self.comp_result_text.insert(tk.END, f"分析模式: {mode}\n")
        if mode == 'flank':
            self.comp_result_text.insert(tk.END, f"Flank分析模式: {self.analyzer.flank_mode}\n")
        else:  # global
            self.comp_result_text.insert(tk.END, f"仅分析GC序列: {'是' if gc_only else '否'}\n")
            
        # 显示结果摘要
        self.comp_result_text.insert(tk.END, "\n互补分析结果:\n")
        self.comp_result_text.insert(tk.END, f"高效率组平均分值: {self.comparison_results['high_efficiency']['mean']:.2f} ± {self.comparison_results['high_efficiency']['std']:.2f}\n")
        self.comp_result_text.insert(tk.END, f"低效率组平均分值: {self.comparison_results['low_efficiency']['mean']:.2f} ± {self.comparison_results['low_efficiency']['std']:.2f}\n")
        self.comp_result_text.insert(tk.END, f"分值差异: {self.comparison_results['difference']['mean_diff']:.2f} ({self.comparison_results['difference']['percent_diff']:.1f}%)\n\n")
        
        # 显示自由能结果
        self.comp_result_text.insert(tk.END, f"高效率组平均自由能: {self.comparison_results['high_efficiency']['energy_mean']:.2f} ± {self.comparison_results['high_efficiency']['energy_std']:.2f} kcal/mol\n")
        self.comp_result_text.insert(tk.END, f"低效率组平均自由能: {self.comparison_results['low_efficiency']['energy_mean']:.2f} ± {self.comparison_results['low_efficiency']['energy_std']:.2f} kcal/mol\n")
        self.comp_result_text.insert(tk.END, f"自由能差异: {self.comparison_results['difference']['energy_mean_diff']:.2f} kcal/mol ({self.comparison_results['difference']['energy_percent_diff']:.1f}%)\n\n")
        
        self.comp_result_text.insert(tk.END, f"图表类型: {self.comp_plot_type_var.get()}\n\n")
        
        self.comp_result_text.insert(tk.END, f"高效率组样本数量: {len(self.comparison_results['high_efficiency']['scores'])}\n")
        self.comp_result_text.insert(tk.END, f"低效率组样本数量: {len(self.comparison_results['low_efficiency']['scores'])}\n\n")
        
        self.comp_result_text.insert(tk.END, "高效率组分值: " + ", ".join([f"{score:.2f}" for score in self.comparison_results['high_efficiency']['scores']]) + "\n\n")
        self.comp_result_text.insert(tk.END, "低效率组分值: " + ", ".join([f"{score:.2f}" for score in self.comparison_results['low_efficiency']['scores']]) + "\n\n")
        
        self.comp_result_text.insert(tk.END, "高效率组自由能: " + ", ".join([f"{energy:.2f}" for energy in self.comparison_results['high_efficiency']['energies']]) + "\n\n")
        self.comp_result_text.insert(tk.END, "低效率组自由能: " + ", ".join([f"{energy:.2f}" for energy in self.comparison_results['low_efficiency']['energies']]) + "\n")
        
        # 显示可视化
        for widget in self.comp_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        # 创建标签页控件用于显示不同的图表
        chart_tabs = ttk.Notebook(self.comp_figure_frame)
        chart_tabs.pack(fill=tk.BOTH, expand=True)
        
        # 创建互补评分图表标签页
        score_tab = ttk.Frame(chart_tabs)
        chart_tabs.add(score_tab, text="互补评分比较")
        
        # 创建自由能图表标签页
        energy_tab = ttk.Frame(chart_tabs)
        chart_tabs.add(energy_tab, text="自由能比较")
        
        # 获取选择的图表类型
        plot_type = self.comp_plot_type_var.get()
        
        # 构建标题，包含分析模式信息
        if mode == 'flank':
            title = f"高低效率序列互补评分比较 (Flank模式: {self.analyzer.flank_mode})"
            energy_title = f"高低效率序列自由能比较 (Flank模式: {self.analyzer.flank_mode})"
        else:  # global
            if gc_only:
                title = "高低效率序列互补评分比较 (全局模式, 仅GC)"
                energy_title = "高低效率序列自由能比较 (全局模式, 仅GC)"
            else:
                title = "高低效率序列互补评分比较 (全局模式)"
                energy_title = "高低效率序列自由能比较 (全局模式)"
        
        # 绘制互补评分比较图表
        plt.figure(figsize=(10, 6))
        self.analyzer.plot_comparison(
            self.comparison_results, 
            title=title,
            plot_type=plot_type
        )
        
        # 将互补评分图表嵌入到第一个标签页
        score_canvas = FigureCanvasTkAgg(plt.gcf(), master=score_tab)
        score_canvas.draw()
        score_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # 绘制自由能比较图表
        plt.figure(figsize=(10, 6))
        self.analyzer.plot_energy_comparison(
            self.comparison_results, 
            title=energy_title,
            plot_type=plot_type
        )
        
        # 将自由能图表嵌入到第二个标签页
        energy_canvas = FigureCanvasTkAgg(plt.gcf(), master=energy_tab)
        energy_canvas.draw()
        energy_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def _clear_compare_inputs(self):
        """清除组比较输入"""
        self.high_eff_path_var.set("")
        self.low_eff_path_var.set("")
        self.output_path_var.set("")
        self.comp_result_text.delete(1.0, tk.END)
        
        # 清除图形区域
        for widget in self.comp_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭所有图形
        plt.close('all')
            
        self.comparison_results = None
        self.status_var.set("就绪")

    def _analyze_a_single_sequence(self):
        """分析单个序列的A碱基susceptibility"""
        sequence = self.a_sequence_text.get(1.0, tk.END).strip()
        if not sequence:
            messagebox.showerror("错误", "请输入DNA序列")
            return
            
        try:
            self.status_var.set("正在分析序列中的A碱基susceptibility...")
            self.root.update_idletasks()
            
            # 保存序列
            self.a_single_sequence = sequence
            
            # 分析序列
            self.a_single_score = self.a_analyzer.analyze_sequence(sequence)
            
            # 显示结果
            self._display_a_single_results(sequence)
            
            self.status_var.set("A碱基susceptibility分析完成")
        except Exception as e:
            messagebox.showerror("错误", str(e))
            self.status_var.set("分析失败")
    
    def _display_a_single_results(self, sequence):
        """显示单序列A碱基susceptibility分析结果"""
        # 清除结果区域
        self.a_result_text.delete(1.0, tk.END)
        
        # 找到序列中所有A碱基的位置
        a_positions = [i for i, base in enumerate(sequence) if base == 'A']
        
        # 显示总评分
        self.a_result_text.insert(tk.END, f"总A碱基susceptibility评分: {self.a_single_score:.2f}\n")
        self.a_result_text.insert(tk.END, f"序列中A碱基数量: {len(a_positions)}\n")
        if len(a_positions) > 0:
            self.a_result_text.insert(tk.END, f"平均每个A碱基分值: {self.a_single_score/len(a_positions):.2f}\n\n")
        
        # 显示详细结果
        if a_positions:
            self.a_result_text.insert(tk.END, "各A碱基susceptibility分析:\n")
            self.a_result_text.insert(tk.END, f"{'位置':<6}{'位置效应':<10}{'嘌呤长度':<10}{'嘌呤效应':<10}{'总分':<10}\n")
            self.a_result_text.insert(tk.END, "-" * 45 + "\n")
            
            for pos in a_positions:
                position_effect = math.log(pos + 1) if pos > 0 else 0
                purine_length = self.a_analyzer._get_continuous_purine_length(sequence, pos)
                purine_effect = purine_length ** 2
                total = position_effect + purine_effect
                
                self.a_result_text.insert(tk.END, f"{pos+1:<6}{position_effect:<10.2f}{purine_length:<10}{purine_effect:<10.2f}{total:<10.2f}\n")
        
        # 显示可视化
        for widget in self.a_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        fig = plt.figure(figsize=(6, 4))
        self.a_analyzer.visualize_sequence_susceptibility(sequence, title=f"序列的A碱基Susceptibility分析")
        
        canvas = FigureCanvasTkAgg(plt.gcf(), master=self.a_figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def _save_a_single_results(self):
        """保存单序列A碱基susceptibility分析结果"""
        if self.a_single_score is None or self.a_single_sequence is None:
            messagebox.showerror("错误", "没有可保存的结果")
            return
            
        filetypes = [
            ("文本文件", "*.txt"),
            ("CSV文件", "*.csv"),
            ("所有文件", "*.*")
        ]
        filename = filedialog.asksaveasfilename(filetypes=filetypes, defaultextension=".txt")
        if not filename:
            return
            
        try:
            # 准备数据
            sequence = self.a_single_sequence
            a_positions = [i for i, base in enumerate(sequence) if base == 'A']
            details = []
            
            for pos in a_positions:
                position_effect = math.log(pos + 1) if pos > 0 else 0
                purine_length = self.a_analyzer._get_continuous_purine_length(sequence, pos)
                purine_effect = purine_length ** 2
                total = position_effect + purine_effect
                
                details.append({
                    'position': pos + 1,
                    'position_effect': position_effect,
                    'purine_length': purine_length,
                    'purine_effect': purine_effect,
                    'score': total
                })
            
            # 保存结果
            ext = os.path.splitext(filename)[1].lower()
            
            if ext == '.csv':
                # 保存为CSV
                df = pd.DataFrame(details)
                df.to_csv(filename, index=False)
            else:
                # 保存为文本
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(f"A碱基Susceptibility分析结果\n")
                    f.write("=" * 50 + "\n\n")
                    
                    f.write(f"总A碱基susceptibility评分: {self.a_single_score:.2f}\n")
                    f.write(f"序列: {sequence}\n")
                    f.write(f"序列长度: {len(sequence)}bp\n")
                    f.write(f"序列中A碱基数量: {len(a_positions)}\n")
                    if len(a_positions) > 0:
                        f.write(f"平均每个A碱基分值: {self.a_single_score/len(a_positions):.2f}\n\n")
                    
                    f.write("各A碱基susceptibility分析:\n")
                    f.write(f"{'位置':<6}{'位置效应':<10}{'嘌呤长度':<10}{'嘌呤效应':<10}{'总分':<10}\n")
                    f.write("-" * 45 + "\n")
                    
                    for detail in details:
                        f.write(f"{detail['position']:<6}{detail['position_effect']:<10.2f}{detail['purine_length']:<10}"
                               f"{detail['purine_effect']:<10.2f}{detail['score']:<10.2f}\n")
            
            # 保存图像
            img_filename = os.path.splitext(filename)[0] + ".png"
            plt.figure(figsize=(8, 6))
            self.a_analyzer.visualize_sequence_susceptibility(sequence, 
                                                         title=f"序列的A碱基Susceptibility分析",
                                                         save_path=img_filename)
            
            self.status_var.set(f"结果已保存至 {filename}")
            messagebox.showinfo("保存成功", f"结果已保存至:\n{filename}\n图像已保存至:\n{img_filename}")
        except Exception as e:
            messagebox.showerror("保存失败", str(e))
    
    def _compare_a_susceptibility(self):
        """比较两组序列的A碱基susceptibility"""
        high_eff_file = self.a_high_eff_path_var.get()
        low_eff_file = self.a_low_eff_path_var.get()
        
        if not high_eff_file or not low_eff_file:
            messagebox.showerror("错误", "请选择高效率和低效率序列文件")
            return
            
        if not os.path.exists(high_eff_file) or not os.path.exists(low_eff_file):
            messagebox.showerror("错误", "文件不存在")
            return
            
        try:
            self.status_var.set("正在读取序列...")
            self.root.update_idletasks()
            
            # 读取序列
            high_eff_seqs = read_sequences_from_file(high_eff_file)
            low_eff_seqs = read_sequences_from_file(low_eff_file)
            
            if not high_eff_seqs or not low_eff_seqs:
                messagebox.showerror("错误", "无法从文件中读取序列")
                self.status_var.set("读取序列失败")
                return
                
            self.status_var.set(f"正在分析 {len(high_eff_seqs)} 条高效率序列和 {len(low_eff_seqs)} 条低效率序列的A碱基susceptibility...")
            self.root.update_idletasks()
            
            # 比较两组
            self.a_comparison_results = self.a_analyzer.compare_groups(high_eff_seqs, low_eff_seqs)
            
            # 显示结果
            self._display_a_comparison_results()
            
            # 保存结果
            output_file = self.a_output_path_var.get()
            if output_file:
                self._save_a_comparison_results(output_file)
                self.status_var.set(f"结果已保存至 {output_file}")
            else:
                self.status_var.set("分析完成")
                
        except Exception as e:
            messagebox.showerror("错误", str(e))
            self.status_var.set("分析失败")
    
    def _display_a_comparison_results(self):
        """显示A碱基susceptibility组比较分析结果"""
        # 清除结果区域
        self.a_comp_result_text.delete(1.0, tk.END)
        
        # 显示结果摘要
        self.a_comp_result_text.insert(tk.END, "A碱基Susceptibility分析结果:\n")
        self.a_comp_result_text.insert(tk.END, f"高效率组平均分值: {self.a_comparison_results['high_efficiency']['mean']:.2f} ± {self.a_comparison_results['high_efficiency']['std']:.2f}\n")
        self.a_comp_result_text.insert(tk.END, f"低效率组平均分值: {self.a_comparison_results['low_efficiency']['mean']:.2f} ± {self.a_comparison_results['low_efficiency']['std']:.2f}\n")
        self.a_comp_result_text.insert(tk.END, f"分值差异: {self.a_comparison_results['difference']['mean_diff']:.2f} ({self.a_comparison_results['difference']['percent_diff']:.1f}%)\n\n")
        
        self.a_comp_result_text.insert(tk.END, f"图表类型: {self.a_plot_type_var.get()}\n\n")
        
        self.a_comp_result_text.insert(tk.END, f"高效率组样本数量: {len(self.a_comparison_results['high_efficiency']['scores'])}\n")
        self.a_comp_result_text.insert(tk.END, f"低效率组样本数量: {len(self.a_comparison_results['low_efficiency']['scores'])}\n\n")
        
        self.a_comp_result_text.insert(tk.END, "高效率组分值: " + ", ".join([f"{score:.2f}" for score in self.a_comparison_results['high_efficiency']['scores']]) + "\n\n")
        self.a_comp_result_text.insert(tk.END, "低效率组分值: " + ", ".join([f"{score:.2f}" for score in self.a_comparison_results['low_efficiency']['scores']]) + "\n")
        
        # 显示可视化
        for widget in self.a_comp_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        # 获取选择的图表类型
        plot_type = self.a_plot_type_var.get()
        
        # 绘制比较图表
        plt.figure(figsize=(10, 6))
        self.a_analyzer.plot_comparison(self.a_comparison_results, 
                                    title="高低效率序列组的A碱基Susceptibility比较",
                                    plot_type=plot_type)
        
        # 嵌入图表
        canvas = FigureCanvasTkAgg(plt.gcf(), master=self.a_comp_figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def _save_a_comparison_results(self, file_path):
        """保存A碱基susceptibility比较结果"""
        if not self.a_comparison_results:
            return
            
        try:
            # 检查文件扩展名
            ext = os.path.splitext(file_path)[1].lower()
            
            if ext == '.json':
                # 保存为JSON
                results_json = {
                    'high_efficiency': {
                        'scores': [float(x) for x in self.a_comparison_results['high_efficiency']['scores']],
                        'mean': float(self.a_comparison_results['high_efficiency']['mean']),
                        'std': float(self.a_comparison_results['high_efficiency']['std'])
                    },
                    'low_efficiency': {
                        'scores': [float(x) for x in self.a_comparison_results['low_efficiency']['scores']],
                        'mean': float(self.a_comparison_results['low_efficiency']['mean']),
                        'std': float(self.a_comparison_results['low_efficiency']['std'])
                    },
                    'difference': {
                        'mean_diff': float(self.a_comparison_results['difference']['mean_diff']),
                        'percent_diff': float(self.a_comparison_results['difference']['percent_diff'])
                    }
                }
                
                with open(file_path, 'w', encoding='utf-8') as f:
                    json.dump(results_json, f, indent=2, ensure_ascii=False)
            
            elif ext == '.csv':
                # 创建CSV报告
                rows = []
                
                # 高效率组
                for i, score in enumerate(self.a_comparison_results['high_efficiency']['scores']):
                    rows.append({
                        'Group': 'High Efficiency Oligo',
                        'Index': i + 1,
                        'A_Susceptibility_Score': score,
                        'Mean_Group_Score': self.a_comparison_results['high_efficiency']['mean'],
                        'Std_Group_Score': self.a_comparison_results['high_efficiency']['std']
                    })
                
                # 低效率组
                for i, score in enumerate(self.a_comparison_results['low_efficiency']['scores']):
                    rows.append({
                        'Group': 'Low Efficiency Oligo',
                        'Index': i + 1,
                        'A_Susceptibility_Score': score,
                        'Mean_Group_Score': self.a_comparison_results['low_efficiency']['mean'],
                        'Std_Group_Score': self.a_comparison_results['low_efficiency']['std']
                    })
                
                # 创建DataFrame并保存
                df = pd.DataFrame(rows)
                df.to_csv(file_path, index=False, encoding='utf-8')
            
            elif ext == '.xlsx':
                # 创建Excel报告
                writer = pd.ExcelWriter(file_path, engine='openpyxl')
                
                # 摘要表
                summary_data = {
                    'Group': ['High Efficiency Oligos', 'Low Efficiency Oligos'],
                    'Mean A Susceptibility Score': [self.a_comparison_results['high_efficiency']['mean'], 
                                                  self.a_comparison_results['low_efficiency']['mean']],
                    'Std Dev': [self.a_comparison_results['high_efficiency']['std'], 
                               self.a_comparison_results['low_efficiency']['std']],
                    'Sample Size': [len(self.a_comparison_results['high_efficiency']['scores']), 
                                   len(self.a_comparison_results['low_efficiency']['scores'])]
                }
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='Summary', index=False)
                
                # 详细分数表
                scores_data = {
                    'High Efficiency Oligo Scores': pd.Series(self.a_comparison_results['high_efficiency']['scores']),
                    'Low Efficiency Oligo Scores': pd.Series(self.a_comparison_results['low_efficiency']['scores'])
                }
                scores_df = pd.DataFrame(scores_data)
                scores_df.to_excel(writer, sheet_name='Scores', index=True)
                
                # 保存
                writer.save()
            
            else:
                # 默认保存为文本
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(f"A碱基Susceptibility分析比较结果\n")
                    f.write("=" * 50 + "\n\n")
                    
                    f.write("高效率组:\n")
                    f.write(f"  样本数量: {len(self.a_comparison_results['high_efficiency']['scores'])}\n")
                    f.write(f"  平均分值: {self.a_comparison_results['high_efficiency']['mean']:.2f} ± {self.a_comparison_results['high_efficiency']['std']:.2f}\n\n")
                    
                    f.write("低效率组:\n")
                    f.write(f"  样本数量: {len(self.a_comparison_results['low_efficiency']['scores'])}\n")
                    f.write(f"  平均分值: {self.a_comparison_results['low_efficiency']['mean']:.2f} ± {self.a_comparison_results['low_efficiency']['std']:.2f}\n\n")
                    
                    f.write("差异分析:\n")
                    f.write(f"  分值差异: {self.a_comparison_results['difference']['mean_diff']:.2f}\n")
                    f.write(f"  差异百分比: {self.a_comparison_results['difference']['percent_diff']:.2f}%\n\n")
                    
                    f.write("详细分值:\n")
                    f.write("  高效率组: " + ", ".join([f"{score:.2f}" for score in self.a_comparison_results['high_efficiency']['scores']]) + "\n")
                    f.write("  低效率组: " + ", ".join([f"{score:.2f}" for score in self.a_comparison_results['low_efficiency']['scores']]) + "\n\n")
            
            # 保存图像
            img_filename = os.path.splitext(file_path)[0] + ".png"
            plt.figure(figsize=(10, 6))
            self.a_analyzer.plot_comparison(self.a_comparison_results, 
                                        title="高低效率序列组的A碱基Susceptibility比较",
                                        plot_type=self.a_plot_type_var.get(),
                                        save_path=img_filename)
            
            messagebox.showinfo("保存成功", f"结果已保存至:\n{file_path}\n图像已保存至:\n{img_filename}")
            
        except Exception as e:
            messagebox.showerror("保存失败", str(e))
    
    def _clear_a_compare_inputs(self):
        """清除A碱基susceptibility比较输入"""
        self.a_high_eff_path_var.set("")
        self.a_low_eff_path_var.set("")
        self.a_output_path_var.set("")
        self.a_comp_result_text.delete(1.0, tk.END)
        
        # 清除图形区域
        for widget in self.a_comp_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭所有图形
        plt.close('all')
            
        self.a_comparison_results = None
        self.status_var.set("就绪")


def save_results_to_file(results: Dict[str, Any], file_path: str) -> None:
    """将分析结果保存到文件"""
    # 检查文件扩展名
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.json':
        # 保存为JSON
        # 转换NumPy数组为列表以便JSON序列化
        results_json = {}
        for group_key, group_data in results.items():
            results_json[group_key] = {}
            for data_key, data_value in group_data.items():
                if data_key in ['scores', 'energies']:
                    results_json[group_key][data_key] = [float(x) for x in data_value]
                elif data_key in ['mean', 'std', 'energy_mean', 'energy_std']:
                    results_json[group_key][data_key] = float(data_value)
                elif data_key == 'details':
                    # 详细结果需要特殊处理
                    details_json = []
                    for detail in data_value:
                        comp_regions, score = detail
                        detail_item = {
                            'score': float(score),
                            'complementary_regions': []
                        }
                        for region in comp_regions:
                            region_copy = {}
                            for k, v in region.items():
                                if isinstance(v, (int, float)):
                                    region_copy[k] = float(v)
                                else:
                                    region_copy[k] = v
                            detail_item['complementary_regions'].append(region_copy)
                        details_json.append(detail_item)
                    results_json[group_key][data_key] = details_json
                else:
                    results_json[group_key][data_key] = data_value
        
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(results_json, f, indent=2, ensure_ascii=False)
    
    elif ext == '.csv':
        # 创建一个简化的CSV报告
        rows = []
        
        # 高效率组
        for i, (score, energy) in enumerate(zip(results['high_efficiency']['scores'], results['high_efficiency']['energies'])):
            rows.append({
                'Group': 'High Efficiency Oligo',
                'Index': i + 1,
                'Score': score,
                'Free_Energy': energy,
                'Mean_Group_Score': results['high_efficiency']['mean'],
                'Std_Group_Score': results['high_efficiency']['std'],
                'Mean_Group_Energy': results['high_efficiency']['energy_mean'],
                'Std_Group_Energy': results['high_efficiency']['energy_std']
            })
        
        # 低效率组
        for i, (score, energy) in enumerate(zip(results['low_efficiency']['scores'], results['low_efficiency']['energies'])):
            rows.append({
                'Group': 'Low Efficiency Oligo',
                'Index': i + 1,
                'Score': score,
                'Free_Energy': energy,
                'Mean_Group_Score': results['low_efficiency']['mean'],
                'Std_Group_Score': results['low_efficiency']['std'],
                'Mean_Group_Energy': results['low_efficiency']['energy_mean'],
                'Std_Group_Energy': results['low_efficiency']['energy_std']
            })
        
        # 创建DataFrame并保存
        df = pd.DataFrame(rows)
        df.to_csv(file_path, index=False, encoding='utf-8')
    
    elif ext == '.xlsx':
        # 创建Excel报告
        writer = pd.ExcelWriter(file_path, engine='openpyxl')
        
        # 摘要表
        summary_data = {
            'Group': ['High Efficiency Oligos', 'Low Efficiency Oligos'],
            'Mean Score': [results['high_efficiency']['mean'], results['low_efficiency']['mean']],
            'Std Dev': [results['high_efficiency']['std'], results['low_efficiency']['std']],
            'Mean Free Energy': [results['high_efficiency']['energy_mean'], results['low_efficiency']['energy_mean']],
            'Energy Std Dev': [results['high_efficiency']['energy_std'], results['low_efficiency']['energy_std']],
            'Sample Size': [len(results['high_efficiency']['scores']), len(results['low_efficiency']['scores'])]
        }
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # 详细分数表
        scores_data = {
            'High Efficiency Oligo Scores': pd.Series(results['high_efficiency']['scores']),
            'Low Efficiency Oligo Scores': pd.Series(results['low_efficiency']['scores']),
            'High Efficiency Oligo Energies': pd.Series(results['high_efficiency']['energies']),
            'Low Efficiency Oligo Energies': pd.Series(results['low_efficiency']['energies'])
        }
        scores_df = pd.DataFrame(scores_data)
        scores_df.to_excel(writer, sheet_name='Scores', index=True)
        
        # 保存
        writer.save()
    
    else:
        # 默认保存为文本报告
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("DNA Self-complementary Structure Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("High Efficiency Oligos Group:\n")
            f.write(f"  Sample Size: {len(results['high_efficiency']['scores'])}\n")
            f.write(f"  Mean Score: {results['high_efficiency']['mean']:.2f} ± {results['high_efficiency']['std']:.2f}\n")
            f.write(f"  Mean Free Energy: {results['high_efficiency']['energy_mean']:.2f} ± {results['high_efficiency']['energy_std']:.2f} kcal/mol\n\n")
            
            f.write("Low Efficiency Oligos Group:\n")
            f.write(f"  Sample Size: {len(results['low_efficiency']['scores'])}\n")
            f.write(f"  Mean Score: {results['low_efficiency']['mean']:.2f} ± {results['low_efficiency']['std']:.2f}\n")
            f.write(f"  Mean Free Energy: {results['low_efficiency']['energy_mean']:.2f} ± {results['low_efficiency']['energy_std']:.2f} kcal/mol\n\n")
            
            f.write("Difference Analysis:\n")
            f.write(f"  Score Difference: {results['difference']['mean_diff']:.2f}\n")
            f.write(f"  Score Percent Difference: {results['difference']['percent_diff']:.2f}%\n")
            f.write(f"  Free Energy Difference: {results['difference']['energy_mean_diff']:.2f} kcal/mol\n")
            f.write(f"  Free Energy Percent Difference: {results['difference']['energy_percent_diff']:.2f}%\n\n")
            
            f.write("Detailed Scores:\n")
            f.write("  High Efficiency Oligos: " + ", ".join([f"{score:.2f}" for score in results['high_efficiency']['scores']]) + "\n")
            f.write("  Low Efficiency Oligos: " + ", ".join([f"{score:.2f}" for score in results['low_efficiency']['scores']]) + "\n\n")
            
            f.write("Detailed Free Energies:\n")
            f.write("  High Efficiency Oligos: " + ", ".join([f"{energy:.2f}" for energy in results['high_efficiency']['energies']]) + "\n")
            f.write("  Low Efficiency Oligos: " + ", ".join([f"{energy:.2f}" for energy in results['low_efficiency']['energies']]) + "\n")


def main():
    try:
        root = tk.Tk()
        app = OligoAnalyzerGUI(root)
        root.mainloop()
    except KeyboardInterrupt:
        print("\n程序被用户中断")
        # 确保清理资源
        plt.close('all')
        if 'root' in locals():
            root.destroy()
    except Exception as e:
        print(f"\n程序发生错误: {e}")
        # 确保清理资源
        plt.close('all')
        if 'root' in locals():
            root.destroy()
    finally:
        print("程序已退出")


if __name__ == "__main__":
    main()