#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DNA寡核苷酸自我互补结构分析工具 - 核心分析模块
"""

import numpy as np
from typing import List, Dict, Tuple, Any, Optional
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio import SeqIO
import os
import json
import pandas as pd
import math


class DNAComplementAnalyzer:
    """DNA序列中自我互补结构分析类"""
    
    def __init__(self, flank_length: int = 20, window_size: int = 10, min_match: int = 4, flank_mode: str = 'both'):
        """
        初始化分析器
        
        参数:
        flank_length: flank区域长度
        window_size: 内部序列滑动窗口大小
        min_match: 最小匹配长度
        flank_mode: flank分析模式，可选值为'both'（两端）、'5p'（仅5'端）或'3p'（仅3'端）
        """
        self.flank_length = flank_length
        self.window_size = window_size
        self.min_match = min_match
        
        # 验证flank_mode参数
        valid_modes = ['both', '5p', '3p']
        if flank_mode not in valid_modes:
            raise ValueError(f"Invalid flank_mode: {flank_mode}. Must be one of {valid_modes}")
        self.flank_mode = flank_mode
    
    def analyze_sequence(self, seq: str) -> Tuple[List[Dict], float]:
        """
        分析单个DNA序列的自我互补结构
        
        参数:
        seq: DNA序列
        
        返回:
        互补区域列表和总评分
        """
        # 标准化序列并检查有效性
        seq = self._validate_sequence(seq)
        seq_length = len(seq)
        
        # 检查序列长度是否足够
        if seq_length < 2 * self.flank_length + self.min_match:
            raise ValueError(f"序列长度({seq_length})太短，无法进行分析。需要至少{2 * self.flank_length + self.min_match}bp")
        
        # 提取5'和3'端的flank区域
        flank_5 = seq[:self.flank_length]
        flank_3 = seq[-self.flank_length:]
        
        # 内部序列
        internal_seq = seq[self.flank_length:seq_length-self.flank_length]
        
        # 分析自我互补
        results = []
        
        # 根据flank_mode决定分析哪一侧的flank
        if self.flank_mode in ['both', '5p']:
            # 检查5'flank与内部序列的互补性
            five_prime_results = self._check_complementarity(
                flank_5, internal_seq, '5p'
            )
            results.extend(five_prime_results)
        
        if self.flank_mode in ['both', '3p']:
            # 检查3'flank与内部序列的互补性
            three_prime_results = self._check_complementarity(
                flank_3, internal_seq, '3p'
            )
            results.extend(three_prime_results)
        
        # 计算总评分
        total_score = sum(result['score'] for result in results)
        
        return results, total_score
    
    def analyze_sequence_global(self, seq: str, gc_only: bool = False) -> Tuple[List[Dict], float]:
        """
        对整个DNA序列进行全局自我互补结构分析
        
        参数:
        seq: DNA序列
        gc_only: 是否只分析GC碱基（如果为True，则过滤掉所有AT碱基）
        
        返回:
        互补区域列表和总评分
        """
        # 标准化序列并检查有效性
        seq = self._validate_sequence(seq)
        
        # 如果需要，过滤掉AT碱基，只保留GC
        if gc_only:
            seq = ''.join([base for base in seq if base in 'GC'])
            if len(seq) < self.min_match:
                raise ValueError(f"过滤AT后的序列长度({len(seq)})太短，无法进行分析。需要至少{self.min_match}bp")
        
        # 使用滑动窗口检查全局互补
        results = self._check_global_complementarity(seq)
        
        # 计算总评分
        total_score = sum(result['score'] for result in results)
        
        return results, total_score
    
    def _check_global_complementarity(self, seq: str) -> List[Dict]:
        """
        检查整个序列内部的互补性（全局分析）
        
        参数:
        seq: 完整序列
        
        返回:
        互补区域列表
        """
        results = []
        seq_length = len(seq)
        
        # 对序列进行滑动窗口分析
        for i in range(seq_length - self.window_size + 1):
            window1 = seq[i:i+self.window_size]
            window1_rev_comp = str(Seq(window1).reverse_complement())
            
            # 查找其余序列中与当前窗口互补的区域
            for j in range(i + self.window_size, seq_length - self.min_match + 1):
                # 尝试不同长度的匹配
                for k in range(self.min_match, min(self.window_size, seq_length - j) + 1):
                    segment2 = seq[j:j+k]
                    
                    # 检查是否存在于反向互补序列中
                    if segment2 in window1_rev_comp:
                        # 计算匹配评分
                        match_score = self._calculate_score(k, segment2, window1)
                        # 计算自由能
                        free_energy = self.calculate_free_energy(segment2, 
                                                                window1[window1_rev_comp.find(segment2):
                                                                      window1_rev_comp.find(segment2)+len(segment2)])
                        
                        # 记录结果
                        results.append({
                            'region1_pos': (i, i+self.window_size),
                            'region2_pos': (j, j+k),
                            'region1_seq': window1,
                            'region2_seq': segment2,
                            'complementary_seq': window1_rev_comp,
                            'match_length': k,
                            'score': match_score,
                            'free_energy': free_energy
                        })
        
        return results
    
    def _validate_sequence(self, seq: str) -> str:
        """验证并标准化DNA序列"""
        seq = seq.upper().strip()
        
        # 检查序列是否只包含有效碱基
        valid_bases = set('ATGCN')
        if not all(base in valid_bases for base in seq):
            invalid_bases = [base for base in seq if base not in valid_bases]
            raise ValueError(f"序列包含无效碱基: {', '.join(invalid_bases)}")
            
        return seq
    
    def _check_complementarity(self, flank: str, internal: str, flank_type: str) -> List[Dict]:
        """
        检查flank区域与内部序列之间的互补性
        
        参数:
        flank: flank区域序列
        internal: 内部序列
        flank_type: flank类型 ('5p' 或 '3p')
        
        返回:
        互补区域列表
        """
        results = []
        
        # 对内部序列进行滑动窗口分析
        for i in range(len(internal) - self.window_size + 1):
            window = internal[i:i+self.window_size]
            window_rev_comp = str(Seq(window).reverse_complement())
            
            # 在flank序列中寻找可能与窗口反向互补匹配的片段
            for j in range(len(flank) - self.min_match + 1):
                # 尝试不同长度的匹配
                for k in range(self.min_match, min(self.window_size, len(flank) - j) + 1):
                    flank_segment = flank[j:j+k]
                    
                    # 检查是否存在于反向互补序列中
                    if flank_segment in window_rev_comp:
                        # 计算匹配评分，传递序列信息以考虑GC/AT权重
                        match_score = self._calculate_score(k, flank_segment, window)
                        # 计算自由能
                        free_energy = self.calculate_free_energy(flank_segment, 
                                                                window[window_rev_comp.find(flank_segment):
                                                                      window_rev_comp.find(flank_segment)+len(flank_segment)])
                        
                        # 记录结果
                        results.append({
                            'flank': flank_type,
                            'flank_pos': (j, j+k),
                            'internal_pos': (i + self.flank_length, i + self.flank_length + self.window_size),
                            'flank_seq': flank_segment,
                            'internal_seq': window,
                            'complementary_seq': window_rev_comp,
                            'match_length': k,
                            'score': match_score,
                            'free_energy': free_energy
                        })
        
        return results
    
    def _calculate_score(self, match_length: int, flank_seq: str = None, internal_seq: str = None) -> float:
        """
        计算互补匹配的评分，考虑匹配长度和GC/AT含量
        
        参数:
        match_length: 匹配长度
        flank_seq: flank区域序列，用于计算GC含量
        internal_seq: 内部序列，用于计算GC含量
        
        返回:
        评分
        """
        # 基本评分公式: 匹配长度的平方 / 窗口大小
        base_score = (match_length ** 2) / self.window_size
        
        # 如果提供了序列，则考虑GC/AT权重
        if flank_seq and internal_seq:
            # 计算GC和AT配对数量
            gc_count = 0
            at_count = 0
            
            # 获取internal_seq的互补序列
            internal_comp = str(Seq(internal_seq).complement())
            
            # 找到flank_seq在internal_comp中的位置
            pos = internal_comp.find(flank_seq)
            if pos >= 0:
                # 提取实际配对的部分
                paired_internal = internal_seq[pos:pos+len(flank_seq)]
                
                for i in range(len(flank_seq)):
                    if i < len(paired_internal):
                        if self._is_complementary(flank_seq[i], paired_internal[i]):
                            if flank_seq[i] in 'GC':
                                gc_count += 1
                            else:
                                at_count += 1
                
                # 应用GC:AT=3:2的权重
                weight_factor = (3 * gc_count + 2 * at_count) / (2 * (gc_count + at_count)) if (gc_count + at_count) > 0 else 1.0
                
                # 调整基础评分
                return base_score * weight_factor
        
        return base_score
    
    def calculate_free_energy(self, seq1: str, seq2: str) -> float:
        """
        使用最近邻热力学模型计算两个序列互补配对的自由能
        
        参数:
        seq1: 第一个序列
        seq2: 第二个序列
        
        返回:
        预估自由能 (kcal/mol)
        """
        # 最近邻热力学参数 (kcal/mol)
        # 这些值来自实验测量的热力学参数
        nn_params = {
            'AA/TT': -1.2, 'AT/TA': -0.9, 'TA/AT': -0.9, 'CA/GT': -1.7,
            'GT/CA': -1.7, 'CT/GA': -1.5, 'GA/CT': -1.5, 'CG/GC': -2.8,
            'GC/CG': -2.8, 'GG/CC': -2.3
        }
        
        # 错配惩罚
        mismatch_penalty = 1.0
        
        # 末端效应修正
        end_penalty = 0.5
        
        # 获取seq2的互补序列
        seq2_comp = str(Seq(seq2).complement())
        min_len = min(len(seq1), len(seq2_comp))
        
        if min_len < 2:
            return 0.0  # 序列太短，无法使用最近邻模型
        
        free_energy = 0.0
        
        # 计算内部二核苷酸步骤的自由能
        for i in range(min_len - 1):
            # 获取当前位置的碱基
            base1_curr = seq1[i]
            base1_next = seq1[i+1]
            base2_curr = seq2_comp[i]
            base2_next = seq2_comp[i+1]
            
            # 检查是否有匹配的碱基对
            match_curr = self._is_complementary(base1_curr, base2_curr)
            match_next = self._is_complementary(base1_next, base2_next)
            
            if match_curr and match_next:
                # 形成二核苷酸步骤
                dinuc = f"{base1_curr}{base1_next}/{base2_curr}{base2_next}"
                
                # 标准化二核苷酸键表示
                if dinuc in nn_params:
                    free_energy += nn_params[dinuc]
                else:
                    # 尝试反向表示
                    rev_dinuc = f"{base1_next}{base1_curr}/{base2_next}{base2_curr}"
                    if rev_dinuc in nn_params:
                        free_energy += nn_params[rev_dinuc]
                    else:
                        # 如果找不到参数，使用平均值
                        free_energy += -1.8
            elif not match_curr and not match_next:
                # 两个连续的错配
                free_energy += mismatch_penalty * 2
            else:
                # 单个错配
                free_energy += mismatch_penalty
        
        # 添加末端效应修正
        free_energy += end_penalty * 2
        
        return free_energy
    
    def _is_complementary(self, base1: str, base2: str) -> bool:
        """
        检查两个碱基是否互补
        
        参数:
        base1: 第一个碱基
        base2: 第二个碱基
        
        返回:
        如果互补则为True，否则为False
        """
        complementary_pairs = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'
        }
        
        return base2 == complementary_pairs.get(base1, 'X')
    
    def batch_analyze(self, sequences: List[str], mode: str = 'flank', gc_only: bool = False) -> List[Tuple[List[Dict], float]]:
        """
        批量分析多个序列
        
        参数:
        sequences: DNA序列列表
        mode: 分析模式，'flank'表示使用flank分析，'global'表示使用全局分析
        gc_only: 当mode为'global'时，是否只分析GC碱基
        
        返回:
        每个序列的分析结果列表
        """
        results = []
        for seq in sequences:
            if mode == 'flank':
                result = self.analyze_sequence(seq)
            elif mode == 'global':
                result = self.analyze_sequence_global(seq, gc_only)
            else:
                raise ValueError(f"无效的分析模式: {mode}. 必须是 'flank' 或 'global'")
            results.append(result)
        return results
    
    def compare_groups(self, 
                      high_efficiency_oligos: List[str], 
                      low_efficiency_oligos: List[str],
                      mode: str = 'flank',
                      gc_only: bool = False) -> Dict[str, Any]:
        """
        比较高效率和低效率oligos的自我互补性
        
        参数:
        high_efficiency_oligos: 高效率oligo序列列表
        low_efficiency_oligos: 低效率oligo序列列表
        mode: 分析模式，'flank'表示使用flank分析，'global'表示使用全局分析
        gc_only: 当mode为'global'时，是否只分析GC碱基
        
        返回:
        比较结果的字典
        """
        # 分析并收集两组oligo的评分
        high_eff_results = self.batch_analyze(high_efficiency_oligos, mode, gc_only)
        low_eff_results = self.batch_analyze(low_efficiency_oligos, mode, gc_only)
        
        high_eff_scores = [result[1] for result in high_eff_results]
        low_eff_scores = [result[1] for result in low_eff_results]
        
        # 计算自由能
        high_eff_energies = []
        low_eff_energies = []
        
        # 计算高效率组的平均自由能
        for result in high_eff_results:
            regions = result[0]
            if regions:
                # 计算所有互补区域的平均自由能
                total_energy = sum(region['free_energy'] for region in regions)
                avg_energy = total_energy / len(regions) if regions else 0
                high_eff_energies.append(avg_energy)
            else:
                high_eff_energies.append(0)  # 如果没有互补区域，自由能为0
        
        # 计算低效率组的平均自由能
        for result in low_eff_results:
            regions = result[0]
            if regions:
                # 计算所有互补区域的平均自由能
                total_energy = sum(region['free_energy'] for region in regions)
                avg_energy = total_energy / len(regions) if regions else 0
                low_eff_energies.append(avg_energy)
            else:
                low_eff_energies.append(0)  # 如果没有互补区域，自由能为0
        
        # 计算统计数据 - 互补评分
        high_mean = np.mean(high_eff_scores)
        high_std = np.std(high_eff_scores)
        low_mean = np.mean(low_eff_scores)
        low_std = np.std(low_eff_scores)
        
        # 计算统计数据 - 自由能
        high_energy_mean = np.mean(high_eff_energies)
        high_energy_std = np.std(high_eff_energies)
        low_energy_mean = np.mean(low_eff_energies)
        low_energy_std = np.std(low_eff_energies)
        
        # 统计分析结果
        comparison_result = {
            'high_efficiency': {
                'scores': high_eff_scores,
                'mean': high_mean,
                'std': high_std,
                'energies': high_eff_energies,
                'energy_mean': high_energy_mean,
                'energy_std': high_energy_std,
                'details': high_eff_results
            },
            'low_efficiency': {
                'scores': low_eff_scores,
                'mean': low_mean,
                'std': low_std,
                'energies': low_eff_energies,
                'energy_mean': low_energy_mean,
                'energy_std': low_energy_std,
                'details': low_eff_results
            },
            'difference': {
                'mean_diff': low_mean - high_mean,
                'percent_diff': ((low_mean - high_mean) / high_mean) * 100 if high_mean > 0 else float('inf'),
                'energy_mean_diff': low_energy_mean - high_energy_mean,
                'energy_percent_diff': ((low_energy_mean - high_energy_mean) / high_energy_mean) * 100 if high_energy_mean != 0 else float('inf')
            }
        }
        
        return comparison_result
    
    def visualize_global_structure(self, seq: str, results: List[Dict], 
                                 gc_only: bool = False,
                                 title: str = 'DNA Sequence Global Self-complementary Structure', 
                                 save_path: Optional[str] = None) -> None:
        """
        可视化DNA序列中的全局自我互补结构
        
        参数:
        seq: DNA序列
        results: 全局自我互补分析结果（来自analyze_sequence_global方法）
        gc_only: 是否只分析了GC碱基
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        """
        if not results:
            print("No complementary structures found, cannot visualize")
            return
        
        # 设置图表
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # 显示原始序列或GC-only序列
        display_seq = seq
        if gc_only:
            display_seq = ''.join([base for base in seq if base in 'GC'])
            title += " (GC Only)"
        
        seq_length = len(display_seq)
        
        # 绘制序列框架
        ax.plot([0, seq_length], [0, 0], 'k-', linewidth=2)
        
        # 绘制互补区域连接
        colors = plt.cm.tab10(np.linspace(0, 1, min(10, len(results))))
        
        for i, result in enumerate(results):
            color = colors[i % len(colors)]
            
            region1_start = result['region1_pos'][0]
            region1_end = result['region1_pos'][1]
            region2_start = result['region2_pos'][0]
            region2_end = result['region2_pos'][1]
            
            # 绘制连接线
            ax.plot([region1_start, region2_start], [-0.5, -1.0], '-', color=color, alpha=0.7)
            ax.plot([region1_end, region2_end], [-0.5, -1.0], '-', color=color, alpha=0.7)
            
            # 高亮显示区域1
            ax.plot([region1_start, region1_end], [-0.5, -0.5], 'o-', color=color, linewidth=2)
            
            # 高亮显示区域2
            ax.plot([region2_start, region2_end], [-1.0, -1.0], 'o-', color=color, linewidth=2)
            
            # 添加信息标签
            midpoint_x = (region1_start + region2_end) / 2
            ax.text(midpoint_x, -1.5, 
                   f"Match {i+1}: {result['match_length']}bp\n"
                   f"Score: {result['score']:.2f}\n"
                   f"ΔG: {result['free_energy']:.1f} kcal/mol",
                   ha='center', va='top', 
                   bbox=dict(boxstyle='round', fc='white', ec=color, alpha=0.7))
        
        # 设置图表属性
        ax.set_xlim(-5, seq_length+5)
        ax.set_ylim(-4, 1)
        ax.set_title(title, fontsize=16)
        ax.set_xlabel('Sequence Position (bp)')
        ax.set_yticks([])
        
        # 添加总结信息
        total_score = sum(result['score'] for result in results)
        plt.figtext(0.5, 0.01, 
                   f"Total Score: {total_score:.2f}   Number of Complementary Regions: {len(results)}", 
                   ha='center', fontsize=12)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        
        # 保存或显示图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()
    
    def plot_comparison(self, comparison_result: Dict[str, Any], 
                        title: str = 'Comparison of High and Low Efficiency Oligos',
                        save_path: Optional[str] = None,
                        plot_type: str = 'boxviolin') -> None:
        """
        可视化比较结果
        
        参数:
        comparison_result: 比较结果字典（来自compare_groups方法）
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        plot_type: 图表类型，可选值为:
                  'boxviolin' - 箱线图和小提琴图 (默认)
                  'histogram' - 直方图
                  'scatter' - 散点图
                  'boxplot' - 仅箱线图
                  'violinplot' - 仅小提琴图
                  'swarmplot' - 蜂群图
                  'stripplot' - 条带图
        """
        # 设置样式
        sns.set_style("whitegrid")
        plt.figure(figsize=(12, 8))
        
        # 准备数据
        high_scores = comparison_result['high_efficiency']['scores']
        low_scores = comparison_result['low_efficiency']['scores']
        
        # 根据plot_type选择不同的图表类型
        if plot_type == 'boxviolin':
            # 创建箱线图
            ax1 = plt.subplot(121)
            sns.boxplot(data=[high_scores, low_scores], ax=ax1)
            ax1.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax1.set_ylabel('Self-complementarity Score')
            ax1.set_title('Box Plot Comparison')
            
            # 创建小提琴图
            ax2 = plt.subplot(122)
            sns.violinplot(data=[high_scores, low_scores], ax=ax2)
            ax2.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax2.set_ylabel('Self-complementarity Score')
            ax2.set_title('Distribution Comparison')
            
        elif plot_type == 'histogram':
            # 创建直方图
            plt.subplot(111)
            plt.hist([high_scores, low_scores], bins=15, alpha=0.7, 
                    label=['High Efficiency', 'Low Efficiency'])
            plt.xlabel('Self-complementarity Score')
            plt.ylabel('Frequency')
            plt.legend()
            plt.title('Score Distribution Histogram')
            
        elif plot_type == 'scatter':
            # 创建散点图
            plt.subplot(111)
            high_x = np.ones(len(high_scores)) * 1
            low_x = np.ones(len(low_scores)) * 2
            plt.scatter(high_x, high_scores, alpha=0.7, label='High Efficiency')
            plt.scatter(low_x, low_scores, alpha=0.7, label='Low Efficiency')
            plt.xticks([1, 2], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Self-complementarity Score')
            plt.legend()
            plt.title('Score Scatter Plot')
            
        elif plot_type == 'boxplot':
            # 仅箱线图
            plt.subplot(111)
            sns.boxplot(data=[high_scores, low_scores])
            plt.xticks([0, 1], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Self-complementarity Score')
            plt.title('Box Plot Comparison')
            
        elif plot_type == 'violinplot':
            # 仅小提琴图
            plt.subplot(111)
            sns.violinplot(data=[high_scores, low_scores])
            plt.xticks([0, 1], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Self-complementarity Score')
            plt.title('Violin Plot Comparison')
            
        elif plot_type == 'swarmplot':
            # 蜂群图
            plt.subplot(111)
            # 创建DataFrame以便使用seaborn的分类图
            data = []
            for score in high_scores:
                data.append({'Group': 'High Efficiency', 'Score': score})
            for score in low_scores:
                data.append({'Group': 'Low Efficiency', 'Score': score})
            df = pd.DataFrame(data)
            
            sns.swarmplot(x='Group', y='Score', data=df)
            plt.ylabel('Self-complementarity Score')
            plt.title('Swarm Plot Comparison')
            
        elif plot_type == 'stripplot':
            # 条带图
            plt.subplot(111)
            # 创建DataFrame以便使用seaborn的分类图
            data = []
            for score in high_scores:
                data.append({'Group': 'High Efficiency', 'Score': score})
            for score in low_scores:
                data.append({'Group': 'Low Efficiency', 'Score': score})
            df = pd.DataFrame(data)
            
            sns.stripplot(x='Group', y='Score', data=df, jitter=True)
            plt.ylabel('Self-complementarity Score')
            plt.title('Strip Plot Comparison')
            
        else:
            raise ValueError(f"Unsupported plot type: {plot_type}. Must be one of: 'boxviolin', 'histogram', 'scatter', 'boxplot', 'violinplot', 'swarmplot', 'stripplot'")
        
        # 添加统计信息
        plt.figtext(0.5, 0.01, 
                   f"High Efficiency Mean: {comparison_result['high_efficiency']['mean']:.2f} ± {comparison_result['high_efficiency']['std']:.2f}\n"
                   f"Low Efficiency Mean: {comparison_result['low_efficiency']['mean']:.2f} ± {comparison_result['low_efficiency']['std']:.2f}\n"
                   f"Difference: {comparison_result['difference']['mean_diff']:.2f} ({comparison_result['difference']['percent_diff']:.1f}%)",
                   ha='center', fontsize=12)
        
        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # 保存图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    def plot_energy_comparison(self, comparison_result: Dict[str, Any], 
                           title: str = 'Free Energy Comparison of High and Low Efficiency Oligos',
                           save_path: Optional[str] = None,
                           plot_type: str = 'boxviolin') -> None:
        """
        可视化自由能比较结果
        
        参数:
        comparison_result: 比较结果字典（来自compare_groups方法）
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        plot_type: 图表类型，可选值为:
                  'boxviolin' - 箱线图和小提琴图 (默认)
                  'histogram' - 直方图
                  'scatter' - 散点图
                  'boxplot' - 仅箱线图
                  'violinplot' - 仅小提琴图
                  'swarmplot' - 蜂群图
                  'stripplot' - 条带图
        """
        # 设置样式
        sns.set_style("whitegrid")
        plt.figure(figsize=(12, 8))
        
        # 准备数据
        high_energies = comparison_result['high_efficiency']['energies']
        low_energies = comparison_result['low_efficiency']['energies']
        
        # 根据plot_type选择不同的图表类型
        if plot_type == 'boxviolin':
            # 创建箱线图
            ax1 = plt.subplot(121)
            sns.boxplot(data=[high_energies, low_energies], ax=ax1)
            ax1.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax1.set_ylabel('Free Energy (kcal/mol)')
            ax1.set_title('Box Plot Comparison')
            
            # 创建小提琴图
            ax2 = plt.subplot(122)
            sns.violinplot(data=[high_energies, low_energies], ax=ax2)
            ax2.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax2.set_ylabel('Free Energy (kcal/mol)')
            ax2.set_title('Distribution Comparison')
            
        elif plot_type == 'histogram':
            # 创建直方图
            plt.subplot(111)
            plt.hist([high_energies, low_energies], bins=15, alpha=0.7, 
                    label=['High Efficiency', 'Low Efficiency'])
            plt.xlabel('Free Energy (kcal/mol)')
            plt.ylabel('Frequency')
            plt.legend()
            plt.title('Energy Distribution Histogram')
            
        elif plot_type == 'scatter':
            # 创建散点图
            plt.subplot(111)
            high_x = np.ones(len(high_energies)) * 1
            low_x = np.ones(len(low_energies)) * 2
            plt.scatter(high_x, high_energies, alpha=0.7, label='High Efficiency')
            plt.scatter(low_x, low_energies, alpha=0.7, label='Low Efficiency')
            plt.xticks([1, 2], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Free Energy (kcal/mol)')
            plt.legend()
            plt.title('Energy Scatter Plot')
            
        elif plot_type == 'boxplot':
            # 仅箱线图
            plt.subplot(111)
            sns.boxplot(data=[high_energies, low_energies])
            plt.xticks([0, 1], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Free Energy (kcal/mol)')
            plt.title('Box Plot Comparison')
            
        elif plot_type == 'violinplot':
            # 仅小提琴图
            plt.subplot(111)
            sns.violinplot(data=[high_energies, low_energies])
            plt.xticks([0, 1], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('Free Energy (kcal/mol)')
            plt.title('Violin Plot Comparison')
            
        elif plot_type == 'swarmplot':
            # 蜂群图
            plt.subplot(111)
            # 创建DataFrame以便使用seaborn的分类图
            data = []
            for energy in high_energies:
                data.append({'Group': 'High Efficiency', 'Energy': energy})
            for energy in low_energies:
                data.append({'Group': 'Low Efficiency', 'Energy': energy})
            df = pd.DataFrame(data)
            
            sns.swarmplot(x='Group', y='Energy', data=df)
            plt.ylabel('Free Energy (kcal/mol)')
            plt.title('Swarm Plot Comparison')
            
        elif plot_type == 'stripplot':
            # 条带图
            plt.subplot(111)
            # 创建DataFrame以便使用seaborn的分类图
            data = []
            for energy in high_energies:
                data.append({'Group': 'High Efficiency', 'Energy': energy})
            for energy in low_energies:
                data.append({'Group': 'Low Efficiency', 'Energy': energy})
            df = pd.DataFrame(data)
            
            sns.stripplot(x='Group', y='Energy', data=df, jitter=True)
            plt.ylabel('Free Energy (kcal/mol)')
            plt.title('Strip Plot Comparison')
            
        else:
            raise ValueError(f"Unsupported plot type: {plot_type}. Must be one of: 'boxviolin', 'histogram', 'scatter', 'boxplot', 'violinplot', 'swarmplot', 'stripplot'")
        
        # 添加统计信息
        plt.figtext(0.5, 0.01, 
                   f"High Efficiency Mean Energy: {comparison_result['high_efficiency']['energy_mean']:.2f} ± {comparison_result['high_efficiency']['energy_std']:.2f} kcal/mol\n"
                   f"Low Efficiency Mean Energy: {comparison_result['low_efficiency']['energy_mean']:.2f} ± {comparison_result['low_efficiency']['energy_std']:.2f} kcal/mol\n"
                   f"Energy Difference: {comparison_result['difference']['energy_mean_diff']:.2f} kcal/mol ({comparison_result['difference']['energy_percent_diff']:.1f}%)",
                   ha='center', fontsize=12)
        
        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # 保存图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    def visualize_structure(self, seq: str, results: List[Dict], 
                           title: str = 'DNA Sequence Self-complementary Structure', 
                           save_path: Optional[str] = None) -> None:
        """
        可视化DNA序列中的自我互补结构
        
        参数:
        seq: DNA序列
        results: 自我互补分析结果（来自analyze_sequence方法）
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        """
        if not results:
            print("No complementary structures found, cannot visualize")
            return
        
        # 设置图表
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # 基本序列信息
        seq_length = len(seq)
        flank_5 = seq[:self.flank_length]
        flank_3 = seq[-self.flank_length:]
        internal = seq[self.flank_length:seq_length-self.flank_length]
        
        # 绘制序列框架
        ax.plot([0, seq_length], [0, 0], 'k-', linewidth=2)
        
        # 根据flank_mode标记flank区域
        if self.flank_mode in ['both', '5p']:
            ax.fill_between([0, self.flank_length], -0.2, 0.2, color='skyblue', alpha=0.5)
            ax.text(self.flank_length/2, 0.3, "5' Flank", ha='center')
        
        if self.flank_mode in ['both', '3p']:
            ax.fill_between([seq_length-self.flank_length, seq_length], -0.2, 0.2, color='skyblue', alpha=0.5)
            ax.text(seq_length-self.flank_length/2, 0.3, "3' Flank", ha='center')
        
        # 添加内部序列标签
        ax.text(self.flank_length + len(internal)/2, 0.3, "Internal Sequence", ha='center')
        
        # 绘制互补区域连接
        colors = plt.cm.tab10(np.linspace(0, 1, min(10, len(results))))
        
        for i, result in enumerate(results):
            color = colors[i % len(colors)]
            
            if result['flank'] == '5p':
                start_pos = result['flank_pos'][0]
                end_pos = result['flank_pos'][1]
                internal_start = result['internal_pos'][0]
                internal_end = result['internal_pos'][1]
            else:  # '3p'
                start_pos = seq_length - self.flank_length + result['flank_pos'][0]
                end_pos = seq_length - self.flank_length + result['flank_pos'][1]
                internal_start = result['internal_pos'][0]
                internal_end = result['internal_pos'][1]
            
            # 绘制连接线
            ax.plot([start_pos, internal_start], [-0.5, -1.0], '-', color=color, alpha=0.7)
            ax.plot([end_pos, internal_end], [-0.5, -1.0], '-', color=color, alpha=0.7)
            
            # 高亮显示flank区域
            ax.plot([start_pos, end_pos], [-0.5, -0.5], 'o-', color=color, linewidth=2)
            
            # 高亮显示内部区域
            ax.plot([internal_start, internal_end], [-1.0, -1.0], 'o-', color=color, linewidth=2)
            
            # 添加信息标签
            midpoint_x = (start_pos + internal_end) / 2
            ax.text(midpoint_x, -1.5, 
                   f"Match {i+1}: {result['match_length']}bp\n"
                   f"Score: {result['score']:.2f}\n"
                   f"ΔG: {result['free_energy']:.1f} kcal/mol",
                   ha='center', va='top', 
                   bbox=dict(boxstyle='round', fc='white', ec=color, alpha=0.7))
        
        # 设置图表属性
        ax.set_xlim(-5, seq_length+5)
        ax.set_ylim(-4, 1)
        ax.set_title(f"{title} (Mode: {self.flank_mode})", fontsize=16)
        ax.set_xlabel('Sequence Position (bp)')
        ax.set_yticks([])
        
        # 添加总结信息
        total_score = sum(result['score'] for result in results)
        plt.figtext(0.5, 0.01, 
                   f"Total Score: {total_score:.2f}   Number of Complementary Regions: {len(results)}", 
                   ha='center', fontsize=12)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        
        # 保存或显示图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()


class ASusceptibilityAnalyzer:
    """分析DNA序列中A碱基的susceptibility"""
    
    def __init__(self):
        """初始化分析器"""
        pass
    
    def analyze_sequence(self, seq: str) -> float:
        """
        分析DNA序列中A碱基的susceptibility
        
        参数:
        seq: DNA序列
        
        返回:
        序列中所有A碱基susceptibility的总和
        """
        # 标准化序列并检查有效性
        seq = self._validate_sequence(seq)
        
        # 找到所有A碱基的位置（以0为起始索引）
        a_positions = [i for i, base in enumerate(seq) if base == 'A']
        
        if not a_positions:
            return 0.0  # 如果序列中没有A碱基，返回0
        
        total_susceptibility = 0.0
        
        # 计算每个A碱基的susceptibility
        for pos in a_positions:
            # 计算位置效应，第一位为0
            position_effect = math.log(pos + 1) if pos > 0 else 0
            
            # 查找当前A碱基所在的连续嘌呤序列(A和G)
            purine_length = self._get_continuous_purine_length(seq, pos)
            
            # 计算嘌呤堆积效应
            purine_effect = purine_length
            
            # 计算该位置A碱基的susceptibility
            a_susceptibility = position_effect * purine_effect
            
            # 累加到总分
            total_susceptibility += a_susceptibility
        
        return total_susceptibility
    
    def _validate_sequence(self, seq: str) -> str:
        """验证并标准化DNA序列"""
        seq = seq.upper().strip()
        
        # 检查序列是否只包含有效碱基
        valid_bases = set('ATGCN')
        if not all(base in valid_bases for base in seq):
            invalid_bases = [base for base in seq if base not in valid_bases]
            raise ValueError(f"序列包含无效碱基: {', '.join(invalid_bases)}")
            
        return seq
    
    def _get_continuous_purine_length(self, seq: str, pos: int) -> int:
        """
        获取A碱基所在的连续嘌呤序列的长度
        
        参数:
        seq: DNA序列
        pos: A碱基的位置
        
        返回:
        包含该A碱基的连续嘌呤序列的长度
        """
        # 确保pos位置是A碱基
        if seq[pos] != 'A':
            raise ValueError(f"位置 {pos} 的碱基不是A，而是 {seq[pos]}")
        
        # 向左寻找连续嘌呤
        left = pos
        while left > 0 and seq[left-1] in 'AG':
            left -= 1
        
        # 向右寻找连续嘌呤
        right = pos
        while right < len(seq) - 1 and seq[right+1] in 'AG':
            right += 1
        
        # 计算连续嘌呤序列的长度
        return right - left + 1
    
    def batch_analyze(self, sequences: List[str]) -> List[float]:
        """
        批量分析多个序列
        
        参数:
        sequences: DNA序列列表
        
        返回:
        每个序列的A碱基susceptibility分数列表
        """
        results = []
        for seq in sequences:
            score = self.analyze_sequence(seq)
            results.append(score)
        return results
    
    def compare_groups(self, 
                      high_efficiency_oligos: List[str], 
                      low_efficiency_oligos: List[str]) -> Dict[str, Any]:
        """
        比较高效率和低效率oligos的A碱基susceptibility
        
        参数:
        high_efficiency_oligos: 高效率oligo序列列表
        low_efficiency_oligos: 低效率oligo序列列表
        
        返回:
        比较结果的字典
        """
        # 分析并收集两组oligo的评分
        high_eff_scores = self.batch_analyze(high_efficiency_oligos)
        low_eff_scores = self.batch_analyze(low_efficiency_oligos)
        
        # 计算统计数据
        high_mean = np.mean(high_eff_scores)
        high_std = np.std(high_eff_scores)
        low_mean = np.mean(low_eff_scores)
        low_std = np.std(low_eff_scores)
        
        # 统计分析结果
        comparison_result = {
            'high_efficiency': {
                'scores': high_eff_scores,
                'mean': high_mean,
                'std': high_std
            },
            'low_efficiency': {
                'scores': low_eff_scores,
                'mean': low_mean,
                'std': low_std
            },
            'difference': {
                'mean_diff': low_mean - high_mean,
                'percent_diff': ((low_mean - high_mean) / high_mean) * 100 if high_mean > 0 else float('inf')
            }
        }
        
        return comparison_result
    
    def plot_comparison(self, comparison_result: Dict[str, Any], 
                        title: str = 'A Susceptibility Comparison',
                        save_path: Optional[str] = None,
                        plot_type: str = 'boxviolin') -> None:
        """
        可视化比较结果
        
        参数:
        comparison_result: 比较结果字典（来自compare_groups方法）
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        plot_type: 图表类型，可选值同DNAComplementAnalyzer的plot_comparison方法
        """
        # 设置样式
        sns.set_style("whitegrid")
        plt.figure(figsize=(12, 8))
        
        # 准备数据
        high_scores = comparison_result['high_efficiency']['scores']
        low_scores = comparison_result['low_efficiency']['scores']
        
        # 根据plot_type选择不同的图表类型
        if plot_type == 'boxviolin':
            # 创建箱线图
            ax1 = plt.subplot(121)
            sns.boxplot(data=[high_scores, low_scores], ax=ax1)
            ax1.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax1.set_ylabel('A Susceptibility Score')
            ax1.set_title('Box Plot Comparison')
            
            # 创建小提琴图
            ax2 = plt.subplot(122)
            sns.violinplot(data=[high_scores, low_scores], ax=ax2)
            ax2.set_xticklabels(['High Efficiency', 'Low Efficiency'])
            ax2.set_ylabel('A Susceptibility Score')
            ax2.set_title('Distribution Comparison')
            
        elif plot_type == 'histogram':
            # 创建直方图
            plt.subplot(111)
            plt.hist([high_scores, low_scores], bins=15, alpha=0.7, 
                    label=['High Efficiency', 'Low Efficiency'])
            plt.xlabel('A Susceptibility Score')
            plt.ylabel('Frequency')
            plt.legend()
            plt.title('Score Distribution Histogram')
            
        elif plot_type == 'scatter':
            # 创建散点图
            plt.subplot(111)
            high_x = np.ones(len(high_scores)) * 1
            low_x = np.ones(len(low_scores)) * 2
            plt.scatter(high_x, high_scores, alpha=0.7, label='High Efficiency')
            plt.scatter(low_x, low_scores, alpha=0.7, label='Low Efficiency')
            plt.xticks([1, 2], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('A Susceptibility Score')
            plt.legend()
            plt.title('Score Scatter Plot')
            
        elif plot_type in ['boxplot', 'violinplot', 'swarmplot', 'stripplot']:
            # 支持其他类型的图表
            plt.subplot(111)
            
            if plot_type == 'boxplot':
                sns.boxplot(data=[high_scores, low_scores])
            elif plot_type == 'violinplot':
                sns.violinplot(data=[high_scores, low_scores])
            else:
                # 创建DataFrame以便使用seaborn的分类图
                data = []
                for score in high_scores:
                    data.append({'Group': 'High Efficiency', 'Score': score})
                for score in low_scores:
                    data.append({'Group': 'Low Efficiency', 'Score': score})
                df = pd.DataFrame(data)
                
                if plot_type == 'swarmplot':
                    sns.swarmplot(x='Group', y='Score', data=df)
                elif plot_type == 'stripplot':
                    sns.stripplot(x='Group', y='Score', data=df, jitter=True)
            
            plt.xticks([0, 1], ['High Efficiency', 'Low Efficiency'])
            plt.ylabel('A Susceptibility Score')
            plt.title(f'{plot_type.capitalize()} Comparison')
            
        else:
            raise ValueError(f"Unsupported plot type: {plot_type}")
        
        # 添加统计信息
        plt.figtext(0.5, 0.01, 
                   f"High Efficiency Mean: {comparison_result['high_efficiency']['mean']:.2f} ± {comparison_result['high_efficiency']['std']:.2f}\n"
                   f"Low Efficiency Mean: {comparison_result['low_efficiency']['mean']:.2f} ± {comparison_result['low_efficiency']['std']:.2f}\n"
                   f"Difference: {comparison_result['difference']['mean_diff']:.2f} ({comparison_result['difference']['percent_diff']:.1f}%)",
                   ha='center', fontsize=12)
        
        plt.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # 保存图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    def visualize_sequence_susceptibility(self, seq: str, 
                                         title: str = 'A Susceptibility Analysis', 
                                         save_path: Optional[str] = None) -> None:
        """
        可视化DNA序列中A碱基的susceptibility
        
        参数:
        seq: DNA序列
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        """
        # 标准化序列
        seq = self._validate_sequence(seq)
        
        # 找到所有A碱基的位置
        a_positions = [i for i, base in enumerate(seq) if base == 'A']
        
        if not a_positions:
            print("序列中没有A碱基，无法可视化")
            return
        
        # 计算每个A碱基的susceptibility
        susceptibilities = []
        purine_lengths = []
        position_effects = []
        purine_effects = []
        
        for pos in a_positions:
            # 计算位置效应
            position_effect = math.log(pos + 1) if pos > 0 else 0
            position_effects.append(position_effect)
            
            # 计算嘌呤堆积效应
            purine_length = self._get_continuous_purine_length(seq, pos)
            purine_lengths.append(purine_length)
            purine_effect = purine_length ** 2
            purine_effects.append(purine_effect)
            
            # 计算总susceptibility
            susceptibilities.append(position_effect + purine_effect)
        
        # 创建图表
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # 序列可视化轨道
        ax.plot([0, len(seq)], [0, 0], 'k-', linewidth=2)
        
        # 高亮显示A碱基位置
        for i, pos in enumerate(a_positions):
            # 计算颜色强度 - 归一化susceptibility
            max_suscept = max(susceptibilities) if susceptibilities else 1
            color_intensity = susceptibilities[i] / max_suscept if max_suscept > 0 else 0
            
            # 高亮显示A碱基
            ax.scatter(pos, 0, s=120, c=[[1, 0, 0, color_intensity]], zorder=3)
            
            # 显示连续嘌呤区域
            purine_length = purine_lengths[i]
            left = max(0, pos - purine_length // 2)
            ax.plot([left, left + purine_length], [0.2, 0.2], 'g-', linewidth=2, alpha=0.7)
            
            # 添加注释标签
            ax.text(pos, 0.4 + (i % 3) * 0.2, 
                   f"A{i+1}:\nPos={pos+1}\nLen={purine_length}\nScore={susceptibilities[i]:.2f}",
                   ha='center', va='bottom',
                   bbox=dict(boxstyle='round', fc='lightyellow', ec='orange', alpha=0.8))
        
        # 在底部添加序列碱基标签
        for i, base in enumerate(seq):
            color = 'red' if base == 'A' else ('green' if base == 'G' else 'blue')
            ax.text(i, -0.3, base, ha='center', va='center', color=color)
        
        # 绘制susceptibility柱状图
        ax2 = ax.twinx()
        ax2.bar(a_positions, susceptibilities, alpha=0.5, width=0.8, color='lightcoral')
        ax2.set_ylabel('Susceptibility Score')
        
        # 设置坐标轴
        ax.set_xlim(-1, len(seq))
        ax.set_ylim(-0.5, 1.5 + (len(a_positions) // 5) * 0.2)
        ax.set_xlabel('Sequence Position (0-indexed)')
        ax.set_title(title, fontsize=16)
        
        # 隐藏y轴刻度
        ax.set_yticks([])
        
        # 添加图例
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        
        legend_elements = [
            Patch(facecolor='red', alpha=0.5, label='A碱基'),
            Line2D([0], [0], color='green', linewidth=2, label='连续嘌呤序列'),
            Patch(facecolor='lightcoral', alpha=0.5, label='Susceptibility分数')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # 添加总结信息
        total_susceptibility = sum(susceptibilities)
        plt.figtext(0.5, 0.01, 
                   f"总A碱基Susceptibility分数: {total_susceptibility:.2f}\n"
                   f"A碱基数量: {len(a_positions)}   平均分值: {total_susceptibility/len(a_positions):.2f} (每个A碱基)",
                   ha='center', fontsize=12)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        
        # 保存或显示图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()


# 如果直接运行此脚本
if __name__ == "__main__":
    # 示例DNA序列
    test_seq = "ACGTACGTACGTATCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGT"
    
    # 创建分析器实例
    analyzer = DNAComplementAnalyzer(flank_length=10, window_size=10, min_match=4)
    
    # 分析序列
    results, score = analyzer.analyze_sequence(test_seq)
    
    # 打印结果
    print(f"总互补评分: {score:.2f}")
    print(f"发现 {len(results)} 个潜在互补区域")
    
    for i, result in enumerate(results):
        print(f"\n互补区域 {i+1}:")
        print(f"Flank: {result['flank']}, 位置: {result['flank_pos']}")
        print(f"内部位置: {result['internal_pos']}")
        print(f"Flank序列: {result['flank_seq']}")
        print(f"内部序列: {result['internal_seq']}")
        print(f"互补序列: {result['complementary_seq']}")
        print(f"匹配长度: {result['match_length']} bp")
        print(f"评分: {result['score']:.2f}")
        print(f"预估自由能: {result['free_energy']:.2f} kcal/mol")
    
    # 可视化结构
    analyzer.visualize_structure(test_seq, results)
    
    print("\n\n")
    print("=" * 50)
    print("测试全局互补分析功能")
    print("=" * 50)
    
    # 使用全局分析模式
    global_results, global_score = analyzer.analyze_sequence_global(test_seq)
    
    # 打印全局分析结果
    print(f"全局分析总评分: {global_score:.2f}")
    print(f"发现 {len(global_results)} 个潜在全局互补区域")
    
    for i, result in enumerate(global_results):
        print(f"\n全局互补区域 {i+1}:")
        print(f"区域1位置: {result['region1_pos']}")
        print(f"区域2位置: {result['region2_pos']}")
        print(f"区域1序列: {result['region1_seq']}")
        print(f"区域2序列: {result['region2_seq']}")
        print(f"互补序列: {result['complementary_seq']}")
        print(f"匹配长度: {result['match_length']} bp")
        print(f"评分: {result['score']:.2f}")
        print(f"预估自由能: {result['free_energy']:.2f} kcal/mol")
    
    # 可视化全局互补结构
    analyzer.visualize_global_structure(test_seq, global_results)
    
    # 测试GC-only模式
    print("\n\n")
    print("=" * 50)
    print("测试GC-only全局互补分析功能")
    print("=" * 50)
    
    gc_results, gc_score = analyzer.analyze_sequence_global(test_seq, gc_only=True)
    
    # 打印GC-only分析结果
    print(f"GC-only分析总评分: {gc_score:.2f}")
    print(f"发现 {len(gc_results)} 个潜在GC互补区域")
    
    for i, result in enumerate(gc_results):
        print(f"\nGC互补区域 {i+1}:")
        print(f"区域1位置: {result['region1_pos']}")
        print(f"区域2位置: {result['region2_pos']}")
        print(f"区域1序列: {result['region1_seq']}")
        print(f"区域2序列: {result['region2_seq']}")
        print(f"互补序列: {result['complementary_seq']}")
        print(f"匹配长度: {result['match_length']} bp")
        print(f"评分: {result['score']:.2f}")
        print(f"预估自由能: {result['free_energy']:.2f} kcal/mol")
    
    # 可视化GC-only互补结构
    analyzer.visualize_global_structure(test_seq, gc_results, gc_only=True)
    
    print("\n\n")
    print("=" * 50)
    print("测试A碱基Susceptibility分析")
    print("=" * 50)
    
    # 创建A碱基susceptibility分析器
    a_analyzer = ASusceptibilityAnalyzer()
    
    # 测试序列
    test_sequences = [
        "CAGGGGAGT",   # 两个A，处于不同长度的嘌呤序列中
        "ATTCCCGAGGAGTCAG",  # 多个A，位于不同位置
        "GATCGATCGAAAAAATCGATCG"  # 连续多个A
    ]
    
    for idx, seq in enumerate(test_sequences):
        print(f"\n测试序列 {idx+1}: {seq}")
        
        # 计算A碱基susceptibility
        score = a_analyzer.analyze_sequence(seq)
        print(f"总A碱基susceptibility分数: {score:.2f}")
        
        # 查找所有A碱基位置
        a_positions = [i for i, base in enumerate(seq) if base == 'A']
        print(f"A碱基位置: {[pos+1 for pos in a_positions]}")
        
        # 计算每个A碱基的susceptibility
        for pos in a_positions:
            position_effect = math.log(pos + 1) if pos > 0 else 0
            purine_length = a_analyzer._get_continuous_purine_length(seq, pos)
            purine_effect = purine_length ** 2
            a_susceptibility = position_effect + purine_effect
            
            print(f"  位置{pos+1}的A: 位置效应={position_effect:.2f}, 嘌呤长度={purine_length}, 嘌呤效应={purine_effect:.2f}, 总分={a_susceptibility:.2f}")
    
    # 可视化第二个序列
    print("\n可视化序列分析结果...")
    a_analyzer.visualize_sequence_susceptibility(test_sequences[1], title=f"序列'{test_sequences[1]}'的A碱基Susceptibility分析")
    
    # 测试两组序列比较
    print("\n\n比较两组序列的A碱基susceptibility...")
    high_efficiency = ["GATCAGTCGTAG", "ACGTGTGTATA", "GGGGAGGGGAC"]
    low_efficiency = ["GAGAGAGAGAGG", "AAAAAAATAAA", "GATCAAAAAGA"]
    
    comparison = a_analyzer.compare_groups(high_efficiency, low_efficiency)
    print(f"高效率组平均分: {comparison['high_efficiency']['mean']:.2f} ± {comparison['high_efficiency']['std']:.2f}")
    print(f"低效率组平均分: {comparison['low_efficiency']['mean']:.2f} ± {comparison['low_efficiency']['std']:.2f}")
    print(f"差异: {comparison['difference']['mean_diff']:.2f} ({comparison['difference']['percent_diff']:.1f}%)")
    
    # 绘制比较结果
    a_analyzer.plot_comparison(comparison, title="高低效率序列的A碱基Susceptibility比较") 