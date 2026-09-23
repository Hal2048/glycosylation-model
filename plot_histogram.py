import ast 
from pathlib import Path  
import matplotlib  
matplotlib.use("Agg")  
import matplotlib.pyplot as plt  


INPUT_FILE = Path(__file__).parent / "distribution_data.txt"
OUTPUT_FILE = Path(__file__).parent / "top10.png"  
TOP_N = 10  


def load_distribution(path):
    """读取 distribution_data.txt，每一行是一条糖型记录。"""
    entries = []  # 用来装所有糖型
    with open(path, "r", encoding="utf-8") as f:  # 打开文件
        for line in f:  # 一行一行读
            line = line.rstrip("\n")  # 去掉行末的换行符
            if not line or line.startswith("#"):  # 空行或注释行跳过
                continue
            parts = line.split("\t")  # 按 Tab 切成若干段
            if len(parts) < 5:  # 字段不够也跳过
                continue
            try:
                tup = ast.literal_eval(parts[2])  # 第 3 段是 "(1,0,2,...)" 这种字符串，还原成元组
            except (ValueError, SyntaxError):
                continue
            entries.append({  # 把需要的信息打包成一个字典存起来
                "rank": int(parts[0]),  # 排名
                "id": int(parts[1]),  # 糖型编号
                "tuple": tup,  # 糖型组成（12 元组）
                "raw": float(parts[3]),  # 原始强度
                "abundance": float(parts[4]),  # 归一化丰度
            })
    return entries  # 把全部记录返回给调用者


def plot_top_n(entries, top_n, save_path):
    """画 Top-N 糖型丰度柱状图，并保存为 PNG。"""
    if not entries:  # 没读到任何数据就退出
        print("没有读到糖型数据。", file=__import__("sys").stderr)
        return

    # 按排名取前 N 个
    top = sorted(entries, key=lambda e: e["rank"])[:top_n]
    ranks = [e["rank"] for e in top]  # 柱子按排名等距排列
    ids = [e["id"] for e in top]  # X 轴刻度要显示糖型 ID
    abundances = [e["abundance"] for e in top]  # Y 轴：丰度

    # 准备一张画布
    fig, ax = plt.subplots(figsize=(max(8, top_n * 0.9), 6))
    # 画柱子：位置按排名，X 轴文字显示 ID
    ax.bar(ranks, abundances, color="#3a78b5", edgecolor="#1f3a5f")

    # 标题和坐标轴文字
    ax.set_xlabel("ID")
    ax.set_ylabel("Normalized Abundance")
    ax.set_title(f"Top {len(top)} Glycoforms")

    # X 轴：刻度位置用排名，下方文字用糖型 ID
    ax.set_xticks(ranks)
    ax.set_xticklabels([str(i) for i in ids], fontsize=10)

    # Y 轴顶端留点空，方便写数值标签
    ax.set_ylim(0, max(abundances) * 1.25 if abundances else 1.0)

    # 在每根柱子顶上写丰度数值
    for x, v in zip(ranks, abundances):
        ax.text(x, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)

    # 保存成 PNG 文件，然后关掉画布
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"图片已保存到：{save_path}")

    # 保存成 PNG 文件，然后关掉画布
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"图片已保存到：{save_path}")


if __name__ == "__main__":
    # 入口：从这里开始运行
    if not INPUT_FILE.exists():  # 输入文件不存在就直接报错退出
        print(f"找不到数据文件：{INPUT_FILE}", file=__import__("sys").stderr)
        raise SystemExit(1)

    entries = load_distribution(INPUT_FILE)  # 第一步：读数据
    print(f"一共读到 {len(entries)} 条糖型记录")
    plot_top_n(entries, TOP_N, OUTPUT_FILE)  # 第二步：画图 + 保存