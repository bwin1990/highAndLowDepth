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
from dna_analyzer import DNAComplementAnalyzer
from typing import List, Dict, Any
import json
import platform
import matplotlib.font_manager as fm

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
        
        # Flank分析模式选择
        ttk.Label(params_frame, text="Flank分析模式:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.flank_mode_var = tk.StringVar(value="both")
        flank_mode_combo = ttk.Combobox(params_frame, textvariable=self.flank_mode_var, width=10)
        flank_mode_combo['values'] = ('both', '5p', '3p')
        flank_mode_combo['state'] = 'readonly'
        flank_mode_combo.grid(row=3, column=1, padx=5, pady=2)
        
        # 添加模式说明
        mode_frame = ttk.Frame(input_frame)
        mode_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(mode_frame, text="模式说明: both=两端flank, 5p=仅5'端flank, 3p=仅3'端flank").pack(anchor=tk.W)
        
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
        
        # Flank分析模式选择
        ttk.Label(params_frame, text="Flank分析模式:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_flank_mode_var = tk.StringVar(value="both")
        comp_flank_mode_combo = ttk.Combobox(params_frame, textvariable=self.comp_flank_mode_var, width=10)
        comp_flank_mode_combo['values'] = ('both', '5p', '3p')
        comp_flank_mode_combo['state'] = 'readonly'
        comp_flank_mode_combo.grid(row=3, column=1, padx=5, pady=2)
        
        # 图表类型选择
        ttk.Label(params_frame, text="图表类型:").grid(row=4, column=0, sticky=tk.W, padx=5, pady=2)
        self.comp_plot_type_var = tk.StringVar(value="boxviolin")
        comp_plot_type_combo = ttk.Combobox(params_frame, textvariable=self.comp_plot_type_var, width=10)
        comp_plot_type_combo['values'] = ('boxviolin', 'histogram', 'scatter', 'boxplot', 'violinplot', 'swarmplot', 'stripplot')
        comp_plot_type_combo['state'] = 'readonly'
        comp_plot_type_combo.grid(row=4, column=1, padx=5, pady=2)
        
        # 添加图表类型说明
        chart_frame = ttk.Frame(input_frame)
        chart_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(chart_frame, text="图表类型: boxviolin=箱线图+小提琴图, histogram=直方图, scatter=散点图").pack(anchor=tk.W)
        
        # 添加模式说明
        mode_frame = ttk.Frame(input_frame)
        mode_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(mode_frame, text="模式说明: both=两端flank, 5p=仅5'端flank, 3p=仅3'端flank").pack(anchor=tk.W)
        
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
            
            # 分析序列
            self.single_results, self.single_score = self.analyzer.analyze_sequence(sequence)
            
            # 显示结果
            self._display_single_results(sequence)
            
            self.status_var.set("分析完成")
        except Exception as e:
            messagebox.showerror("错误", str(e))
            self.status_var.set("分析失败")
            
    def _display_single_results(self, sequence):
        """显示单序列分析结果"""
        # 清除结果区域
        self.result_text.delete(1.0, tk.END)
        
        # 显示总评分
        self.result_text.insert(tk.END, f"总互补评分: {self.single_score:.2f}\n")
        self.result_text.insert(tk.END, f"发现 {len(self.single_results)} 个潜在互补区域\n")
        self.result_text.insert(tk.END, f"Flank分析模式: {self.analyzer.flank_mode}\n\n")
        
        # 显示详细结果
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
            
        # 显示可视化
        for widget in self.figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        fig = plt.figure(figsize=(6, 4))
        self.analyzer.visualize_structure(sequence, self.single_results)
        
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
                
                f.write(f"总互补评分: {self.single_score:.2f}\n")
                f.write(f"发现 {len(self.single_results)} 个潜在互补区域\n\n")
                
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
                    
            # 保存图像
            img_filename = os.path.splitext(filename)[0] + ".png"
            fig = plt.figure(figsize=(8, 6))
            sequence = self.sequence_text.get(1.0, tk.END).strip()
            self.analyzer.visualize_structure(sequence, self.single_results)
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
            
            # 比较两组
            self.comparison_results = self.analyzer.compare_groups(high_eff_seqs, low_eff_seqs)
            
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
        
        # 显示结果摘要
        self.comp_result_text.insert(tk.END, "Analysis Results:\n")
        self.comp_result_text.insert(tk.END, f"High Efficiency Group Mean Score: {self.comparison_results['high_efficiency']['mean']:.2f} ± {self.comparison_results['high_efficiency']['std']:.2f}\n")
        self.comp_result_text.insert(tk.END, f"Low Efficiency Group Mean Score: {self.comparison_results['low_efficiency']['mean']:.2f} ± {self.comparison_results['low_efficiency']['std']:.2f}\n")
        self.comp_result_text.insert(tk.END, f"Score Difference: {self.comparison_results['difference']['mean_diff']:.2f} ({self.comparison_results['difference']['percent_diff']:.1f}%)\n")
        self.comp_result_text.insert(tk.END, f"Flank Analysis Mode: {self.analyzer.flank_mode}\n")
        self.comp_result_text.insert(tk.END, f"Chart Type: {self.comp_plot_type_var.get()}\n\n")
        
        self.comp_result_text.insert(tk.END, f"High Efficiency Group Sample Size: {len(self.comparison_results['high_efficiency']['scores'])}\n")
        self.comp_result_text.insert(tk.END, f"Low Efficiency Group Sample Size: {len(self.comparison_results['low_efficiency']['scores'])}\n\n")
        
        self.comp_result_text.insert(tk.END, "High Efficiency Group Scores: " + ", ".join([f"{score:.2f}" for score in self.comparison_results['high_efficiency']['scores']]) + "\n\n")
        self.comp_result_text.insert(tk.END, "Low Efficiency Group Scores: " + ", ".join([f"{score:.2f}" for score in self.comparison_results['low_efficiency']['scores']]) + "\n")
        
        # 显示可视化
        for widget in self.comp_figure_frame.winfo_children():
            widget.destroy()
            
        # 关闭之前的图形
        plt.close('all')
        
        # 创建图形并调用plot_comparison
        plt.figure(figsize=(10, 6))
        
        # 获取选择的图表类型
        plot_type = self.comp_plot_type_var.get()
        
        # 调用比较分析方法，传入图表类型
        self.analyzer.plot_comparison(self.comparison_results, plot_type=plot_type)
        
        # 将图形嵌入到Tkinter界面
        canvas = FigureCanvasTkAgg(plt.gcf(), master=self.comp_figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
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
                if data_key == 'scores':
                    results_json[group_key][data_key] = [float(x) for x in data_value]
                elif data_key in ['mean', 'std']:
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
        for i, score in enumerate(results['high_efficiency']['scores']):
            rows.append({
                'Group': 'High Efficiency Oligo',
                'Index': i + 1,
                'Score': score,
                'Mean_Group_Score': results['high_efficiency']['mean'],
                'Std_Group_Score': results['high_efficiency']['std']
            })
        
        # 低效率组
        for i, score in enumerate(results['low_efficiency']['scores']):
            rows.append({
                'Group': 'Low Efficiency Oligo',
                'Index': i + 1,
                'Score': score,
                'Mean_Group_Score': results['low_efficiency']['mean'],
                'Std_Group_Score': results['low_efficiency']['std']
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
            'Sample Size': [len(results['high_efficiency']['scores']), len(results['low_efficiency']['scores'])]
        }
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # 详细分数表
        scores_data = {
            'High Efficiency Oligo Scores': pd.Series(results['high_efficiency']['scores']),
            'Low Efficiency Oligo Scores': pd.Series(results['low_efficiency']['scores'])
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
            f.write(f"  Mean Score: {results['high_efficiency']['mean']:.2f} ± {results['high_efficiency']['std']:.2f}\n\n")
            
            f.write("Low Efficiency Oligos Group:\n")
            f.write(f"  Sample Size: {len(results['low_efficiency']['scores'])}\n")
            f.write(f"  Mean Score: {results['low_efficiency']['mean']:.2f} ± {results['low_efficiency']['std']:.2f}\n\n")
            
            f.write("Difference Analysis:\n")
            f.write(f"  Score Difference: {results['difference']['mean_diff']:.2f}\n")
            f.write(f"  Percent Difference: {results['difference']['percent_diff']:.2f}%\n\n")
            
            f.write("Detailed Scores:\n")
            f.write("  High Efficiency Oligos: " + ", ".join([f"{score:.2f}" for score in results['high_efficiency']['scores']]) + "\n")
            f.write("  Low Efficiency Oligos: " + ", ".join([f"{score:.2f}" for score in results['low_efficiency']['scores']]) + "\n")


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