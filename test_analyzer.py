#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DNA自我互补结构分析工具的测试脚本
"""

import os
import sys
from dna_analyzer import DNAComplementAnalyzer


def test_single_sequence():
    """测试单序列分析功能"""
    print("测试单序列分析...")
    
    # 创建分析器实例
    analyzer = DNAComplementAnalyzer(flank_length=10, window_size=10, min_match=4)
    
    # 测试序列 (此序列故意设计为具有自我互补结构)
    test_seq = "ACGTACGTACGCCCGCATGGCATGCGGGACGTACGTACGT"
    
    # 分析序列
    results, score = analyzer.analyze_sequence(test_seq)
    
    # 验证结果
    print(f"总互补评分: {score:.2f}")
    print(f"发现 {len(results)} 个潜在互补区域")
    
    if len(results) > 0:
        print("单序列分析测试通过 ✓")
    else:
        print("单序列分析测试失败 ✗")
        print("预期找到至少一个互补区域，但未找到")
    
    return len(results) > 0


def test_group_comparison():
    """测试组间比较功能"""
    print("\n测试组间比较...")
    
    # 创建具有不同互补特性的序列组
    high_eff_seqs = [
        "ACGTACGTACGTATCGATCGATCGATCGATCGACGTACGTACGT",
        "ACGTACGTACGTGCATCGATAGCTAGCTAGCTACGTACGTACGT"
    ]
    
    low_eff_seqs = [
        "ACGTACGTACGTGCATCGATGCATCGATACGTACGTACGT",  # 有意设计为具有更多互补
        "ACGTACGTACGTCGAATTCCGAATTCGACGTACGTACGT"     # 有意设计为具有更多互补
    ]
    
    # 创建分析器
    analyzer = DNAComplementAnalyzer(flank_length=10, window_size=10, min_match=4)
    
    # 比较两组
    comparison_results = analyzer.compare_groups(high_eff_seqs, low_eff_seqs)
    
    # 检查结果是否符合预期 (低效率组应当有更高的互补性评分)
    high_mean = comparison_results['high_efficiency']['mean']
    low_mean = comparison_results['low_efficiency']['mean']
    
    print(f"高效率组平均评分: {high_mean:.2f}")
    print(f"低效率组平均评分: {low_mean:.2f}")
    print(f"差异: {comparison_results['difference']['mean_diff']:.2f}")
    
    # 我们预期低效率组分数更高 (低效率组设计为有更多互补)
    expected_result = low_mean > high_mean
    
    if expected_result:
        print("组间比较测试通过 ✓")
    else:
        print("组间比较测试失败 ✗")
        print("预期低效率组有更高的互补评分，但结果相反")
    
    return expected_result


def test_visualization():
    """测试可视化功能"""
    print("\n测试可视化功能...")
    
    # 创建分析器实例
    analyzer = DNAComplementAnalyzer(flank_length=10, window_size=10, min_match=4)
    
    # 测试序列
    test_seq = "ACGTACGTACGCCCGCATGGCATGCGGGACGTACGTACGT"
    
    # 分析序列
    results, _ = analyzer.analyze_sequence(test_seq)
    
    # 尝试保存图像
    test_img = "test_visual.png"
    try:
        analyzer.visualize_structure(test_seq, results, save_path=test_img)
        success = os.path.exists(test_img)
        
        if success:
            print(f"图像已成功保存到 {test_img}")
            print("可视化测试通过 ✓")
        else:
            print("图像保存失败")
            print("可视化测试失败 ✗")
        
        # 清理
        if os.path.exists(test_img):
            os.remove(test_img)
        
        return success
    
    except Exception as e:
        print(f"可视化测试出现错误: {str(e)}")
        print("可视化测试失败 ✗")
        return False


def run_all_tests():
    """运行所有测试"""
    print("=" * 50)
    print("开始测试DNA自我互补结构分析工具")
    print("=" * 50)
    
    # 运行测试
    test1 = test_single_sequence()
    test2 = test_group_comparison()
    test3 = test_visualization()
    
    # 汇总结果
    print("\n" + "=" * 50)
    print("测试结果汇总:")
    print(f"单序列分析: {'通过 ✓' if test1 else '失败 ✗'}")
    print(f"组间比较: {'通过 ✓' if test2 else '失败 ✗'}")
    print(f"可视化功能: {'通过 ✓' if test3 else '失败 ✗'}")
    
    total_tests = 3
    passed_tests = sum([test1, test2, test3])
    
    print(f"\n总测试数: {total_tests}")
    print(f"通过测试: {passed_tests}")
    print(f"失败测试: {total_tests - passed_tests}")
    
    if passed_tests == total_tests:
        print("\n所有测试通过! 工具功能正常 ✓")
        return 0
    else:
        print("\n部分测试失败! 请检查错误 ✗")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests()) 