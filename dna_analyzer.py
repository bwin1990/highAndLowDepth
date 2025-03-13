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


class DNAComplementAnalyzer:
    """DNA序列中自我互补结构分析类"""
    
    def __init__(self, flank_length: int = 20, window_size: int = 10, min_match: int = 4):
        """
        初始化分析器
        
        参数:
        flank_length: flank区域长度
        window_size: 内部序列滑动窗口大小
        min_match: 最小匹配长度
        """
        self.flank_length = flank_length
        self.window_size = window_size
        self.min_match = min_match
    
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
        
        # 检查5'flank与内部序列的互补性
        five_prime_results = self._check_complementarity(
            flank_5, internal_seq, '5p'
        )
        results.extend(five_prime_results)
        
        # 检查3'flank与内部序列的互补性
        three_prime_results = self._check_complementarity(
            flank_3, internal_seq, '3p'
        )
        results.extend(three_prime_results)
        
        # 计算总评分
        total_score = sum(result['score'] for result in results)
        
        return results, total_score
    
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
                        # 计算匹配评分
                        match_score = self._calculate_score(k)
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
    
    def _calculate_score(self, match_length: int) -> float:
        """
        计算互补匹配的评分
        
        参数:
        match_length: 匹配长度
        
        返回:
        评分
        """
        # 基本评分公式: 匹配长度的平方 / 窗口大小
        # 这样设计使得较长的匹配获得更高的权重
        score = (match_length ** 2) / self.window_size
        return score
    
    def calculate_free_energy(self, seq1: str, seq2: str) -> float:
        """
        简化版计算两个序列互补配对的自由能
        
        参数:
        seq1: 第一个序列
        seq2: 第二个序列
        
        返回:
        预估自由能 (kcal/mol)
        """
        # 计算GC和AT配对数量
        gc_count = 0
        at_count = 0
        
        seq2_comp = str(Seq(seq2).complement())
        
        for i in range(min(len(seq1), len(seq2_comp))):
            if (seq1[i] == 'G' and seq2_comp[i] == 'G') or (seq1[i] == 'C' and seq2_comp[i] == 'C'):
                gc_count += 1
            elif (seq1[i] == 'A' and seq2_comp[i] == 'A') or (seq1[i] == 'T' and seq2_comp[i] == 'T'):
                at_count += 1
        
        # 简化的自由能计算 (kcal/mol)
        # GC对贡献约-3 kcal/mol，AT对贡献约-2 kcal/mol
        free_energy = -3 * gc_count - 2 * at_count
        
        return free_energy
    
    def batch_analyze(self, sequences: List[str]) -> List[Tuple[List[Dict], float]]:
        """
        批量分析多个序列
        
        参数:
        sequences: DNA序列列表
        
        返回:
        每个序列的分析结果列表
        """
        results = []
        for seq in sequences:
            result = self.analyze_sequence(seq)
            results.append(result)
        return results
    
    def compare_groups(self, 
                      high_efficiency_oligos: List[str], 
                      low_efficiency_oligos: List[str]) -> Dict[str, Any]:
        """
        比较高效率和低效率oligos的自我互补性
        
        参数:
        high_efficiency_oligos: 高效率oligo序列列表
        low_efficiency_oligos: 低效率oligo序列列表
        
        返回:
        比较结果的字典
        """
        # 分析并收集两组oligo的评分
        high_eff_results = self.batch_analyze(high_efficiency_oligos)
        low_eff_results = self.batch_analyze(low_efficiency_oligos)
        
        high_eff_scores = [result[1] for result in high_eff_results]
        low_eff_scores = [result[1] for result in low_eff_results]
        
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
                'std': high_std,
                'details': high_eff_results
            },
            'low_efficiency': {
                'scores': low_eff_scores,
                'mean': low_mean,
                'std': low_std,
                'details': low_eff_results
            },
            'difference': {
                'mean_diff': low_mean - high_mean,
                'percent_diff': ((low_mean - high_mean) / high_mean) * 100 if high_mean > 0 else float('inf')
            }
        }
        
        return comparison_result
    
    def plot_comparison(self, comparison_result: Dict[str, Any], 
                        title: str = 'Comparison of High and Low Efficiency Oligos',
                        save_path: Optional[str] = None) -> None:
        """
        可视化比较结果
        
        参数:
        comparison_result: 比较结果字典（来自compare_groups方法）
        title: 图表标题
        save_path: 图表保存路径（如果不为None）
        """
        # 设置样式
        sns.set_style("whitegrid")
        plt.figure(figsize=(12, 8))
        
        # 准备数据
        high_scores = comparison_result['high_efficiency']['scores']
        low_scores = comparison_result['low_efficiency']['scores']
        
        # 创建箱线图
        ax = plt.subplot(121)
        sns.boxplot(data=[high_scores, low_scores], ax=ax)
        ax.set_xticklabels(['High Efficiency', 'Low Efficiency'])
        ax.set_ylabel('Self-complementarity Score')
        ax.set_title('Box Plot Comparison')
        
        # 创建小提琴图
        ax = plt.subplot(122)
        sns.violinplot(data=[high_scores, low_scores], ax=ax)
        ax.set_xticklabels(['High Efficiency', 'Low Efficiency'])
        ax.set_ylabel('Self-complementarity Score')
        ax.set_title('Distribution Comparison')
        
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
        
        # 标记flank区域
        ax.fill_between([0, self.flank_length], -0.2, 0.2, color='skyblue', alpha=0.5)
        ax.fill_between([seq_length-self.flank_length, seq_length], -0.2, 0.2, color='skyblue', alpha=0.5)
        
        # 添加文字标签
        ax.text(self.flank_length/2, 0.3, "5' Flank", ha='center')
        ax.text(seq_length-self.flank_length/2, 0.3, "3' Flank", ha='center')
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