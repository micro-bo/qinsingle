#!/usr/bin/env python
# coding: utf-8
## -------------------------------------------------------------------- ##
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

import json
import random
import time
import hashlib
import requests
from pandas.io.json import json_normalize

## -------------------------------------------------------------------- ##
# argparse初始化
parser = argparse.ArgumentParser()
# argparse添加参数
parser.add_argument('--h5ad_list', required = True, dest = "h5ad_list",nargs="+", type=str,help = "要处理的数据库h5ad格式数据集路径,split by ' '")
parser.add_argument('--kingdom', required = True, dest = "kingdom",choices = ['animals', 'plants'],default = 'animals', help = "界，动物（animals）或植物（plants），决定了metadata与standard_celltype去哪里调取接口")
parser.add_argument('--outdir', required = True, dest = "outdir",help = "指定结果文件的输出路径")
args = parser.parse_args()
## -------------------------------------------------------------------- ##
h5ad_list = args.h5ad_list
kingdom = args.kingdom
outdir = args.outdir

## -------------------------------------------------------------------- ##

# 预设定义

# store_raw_count: T,N,是否存储 raw count 数据
summary_stat_df = pd.DataFrame(
    columns=[
        "rds_id",
        "species",
        "store_raw_count",
        "cell_count",
        "gene_count",
        "standard_symbol_count",
        "standard_symbol_propotion",
        "modified_wrong_format_symbol_count",
        "modified_alias_symbol_count",
        "unmodify_alias_symbols",
        "final_duplicate_symbol_count",
    ]
)


# 错误格式转换
gene_wrong_format_conversion_species = {
    "Homo sapiens": "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/human_allgene_inrds_new.xls"
}


# 进行别名原名转换的物种
gene_alias_conversion_species = {
    "Homo sapiens": "/SGRNJ03/pipline_test/huangxinpeng/01.process_groups/h5ad_format_change/HUGO_alias_symbols_drop_duplicates_20230315.csv"
}

# 检测白名单的物种
check_white_list_species = {
    "Homo sapiens": [
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/homo_sapiens_hugo_20230313.xls",
    ],
    "Mus musculus": [
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/mus_musculus_uniq.xls",
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/mus_musculus_nogenename.xls",
    ],
    "Macaca fascicularis": [
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/macaca_fascicularis_uniq.xls",
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/macaca_fascicularis_nogenename.xls",
    ],
    "Chlorocebus sabaeus": [
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/chlorocebus_sabaeus_uniq.xls",
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/chlorocebus_sabaeus_nogenename.xls",
    ],
    "Callithrix jacchus": [
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/callithrix_jacchus_uniq.xls",
        "/SGRNJ01/Public/Script/datascience/pipeline_public/CheckGene/data/tochenss/callithrix_jacchus_nogenename.xls",
    ],
    'Arabidopsis thaliana':[
        "/SGRNJ01/Public/Script/datascience/pipeline_public/convert_public_h5ad/scripts/gene_whitelist/arabidopsis_thaliana_genlist.xls"
    ],
}

annot_to_abbr_ref = "/SGRNJ03/pipline_test/huangxinpeng/01.process_groups/h5ad_format_change/celltype_20230314.csv"


def get_celltype_df(kingdom):
    if kingdom == "animals":
        url = "http://api-bi.internal.singleronbio.com/api/database-api/dock/getCellTypeAnimal"
        appKey = "BoZk3yn0ah"
        appId = "animalSingleron"
        appSecret = "sT7ADhkq7sMXOHlfQ0lSvlfQDAWSmzFzTlmttE9iG"

        timeStamp = str(int(time.time()))
        sign_temp = (
            "&appKey="
            + appKey
            + "&appId="
            + appId
            + "&timeStamp="
            + timeStamp
            + "&appSecret="
            + appSecret
        )
        sign = hashlib.md5(sign_temp.encode()).hexdigest().upper()

        params = {
            "appKey": appKey,
            "appId": appId,
            "timeStamp": timeStamp,
            "sign": sign,
        }

        response = requests.get(url, params=params)
        ## ------------------------------------------------------------------------- ##
        # 返回的结果内容解析为json格式
        res = json.loads(response.text)
        celltype_info = res["data"]
        celltype_info_df = pd.json_normalize(celltype_info)
        # 将所有列转换为category类型
        celltype_info_df = celltype_info_df.astype("category")
        return celltype_info_df
        ## ------------------------------------------------------------------------- ##
    elif kingdom == "plants":
        url = "https://api.alpha.plant.singleronbio.com/api/database-api/dock/getCellTypePlant"
        appKey = "udfhuenieh"
        appId = "plantSingleron"
        appSecret = "NNUgZZncQqLt2A54KaLgNCG3kN9FUSHefGDGY9JFY"

        timeStamp = str(int(time.time()))

        sign_temp = (
            "&appKey="
            + appKey
            + "&appId="
            + appId
            + "&timeStamp="
            + timeStamp
            + "&appSecret="
            + appSecret
        )
        sign = hashlib.md5(sign_temp.encode()).hexdigest().upper()

        params = {
            "appKey": appKey,
            "appId": appId,
            "timeStamp": timeStamp,
            "sign": sign,
        }
        response = requests.get(url, params=params)
        ## ------------------------------------------------------------------------- ##
        # 返回的结果内容解析为json格式
        res = json.loads(response.text)
        celltype_info = res["data"]
        celltype_info_df = pd.json_normalize(celltype_info)
        # 将所有列转换为category类型
        celltype_info_df = celltype_info_df.astype("category")
        return celltype_info_df
    else:
        pass
    ## ------------------------------------------------------------------------- ##



## -------------------------------------------------------------------- ##
# 全称简称转换字典
def convert_abbr():
    convert_name = get_celltype_df(kingdom)
    convert_name = convert_name[["cellTypeStandard", "abbreviation"]]
    convert_name = convert_name[convert_name["abbreviation"].notna()].reset_index(
        drop=True
    )

    convert_abbr_dict, convert_reverse_dict = {}, {}
    for name in range(len(convert_name.index)):
        convert_abbr_dict[convert_name["cellTypeStandard"][name]] = convert_name[
            "abbreviation"
        ][name]
        convert_reverse_dict[convert_name["abbreviation"][name]] = convert_name[
            "cellTypeStandard"
        ][name]

    # cluster_standard 中有而不在ref中的全称
    annot_diff = list(
        set(adata.obs["cluster_standard"].cat.categories).difference(
            set(convert_name["cellTypeStandard"])
        )
    )
    if len(annot_diff) > 0:
        #    for name2 in range(len(annot_diff)):
        #        convert_abbr_dict[str(annot_diff[name2])] = str(annot_diff[name2]).title().replace(" ", "")
        return sys.exit(
            "error: h5ad中存在未在DocManage中出现的新细胞类型:{}，暂未定义其细胞类型简称！".format(annot_diff)
        )
    else:
        return convert_abbr_dict, convert_reverse_dict

## -------------------------------------------------------------------- ##
# 循环输入的数据集列表
for a in range(0, len(h5ad_list), 1):
    adata = sc.read(h5ad_list[a])

    species = adata.obs["Species"].cat.categories[0]
    rds_id = os.path.basename(h5ad_list[a]).rsplit("_", 1)[0]
    print("——————正在处理: ", rds_id, " ——————")
    cell_count = adata.n_obs
    # 统计基因数
    ori_gene_name = list(adata.var.index)  # 原始基因名
    gene_count = len(ori_gene_name)  # 基因数
    modified_symbols = ori_gene_name  # 不一定修改
    ## -------------------------------------------------------------------- ##
    # 如果是需要基因名格式转换的物种，则读取对应的转换表，并转成字典
    if species in gene_wrong_format_conversion_species:
        wrong_format_symbols_df = pd.read_table(
            gene_wrong_format_conversion_species[species]
        )

        # get the column names as a list
        cols = wrong_format_symbols_df.columns.tolist()
        # rename the first and second columns
        cols[0] = "ori_symbols"
        cols[1] = "new_symbols"
        wrong_format_symbols_df.columns = cols

        wrong_format_symbols_dict = wrong_format_symbols_df.set_index("ori_symbols")[
            "new_symbols"
        ].to_dict()

        print("数据集包含基因数：", gene_count)
        # 仅需要格式转换的物种计算
        need_modify_wrong_format_symbols = set(ori_gene_name) & set(
            wrong_format_symbols_df["ori_symbols"]
        )
        modify_wrong_format_symbol_count = len(need_modify_wrong_format_symbols)
        print(
            "原始基因名中有",
            modify_wrong_format_symbol_count,
            "个基因需要进行错误格式修改:",
            need_modify_wrong_format_symbols,
        )
    else:
        modify_wrong_format_symbol_count = -1
    ## -------------------------------------------------------------------- ##
    # 如果是需要别名转换的物种，则读取对应的转换表，并转成字典
    if species in gene_alias_conversion_species:
        HUGO_alias_symbols_df = pd.read_csv(gene_alias_conversion_species[species])

        # get the column names as a list
        cols = HUGO_alias_symbols_df.columns.tolist()
        # rename the first and second columns
        cols[0] = "ori_symbols"
        cols[1] = "new_symbols"
        HUGO_alias_symbols_df.columns = cols

        HUGO_alias_symbols_dict = HUGO_alias_symbols_df.set_index("ori_symbols")[
            "new_symbols"
        ].to_dict()

        # 仅需要别名转换的物种计算
        need_modify_alias_symbols = set(modified_symbols) & set(
            HUGO_alias_symbols_df["ori_symbols"]
        )

        modify_alias_symbols = len(need_modify_alias_symbols)
        print(
            "（经历错误格式转换后，）有",
            modify_alias_symbols,
            "个基因需要进行别名修改:",
            ", ".join(need_modify_alias_symbols)[:1000],
            "（截取前1000个字符，最后一个基因名不一定完整）",
        )

        if len(need_modify_alias_symbols) > 0:
            # 先尝试转换看有几个转换后重复的基因
            try_modified_symbols = [
                HUGO_alias_symbols_dict.get(gene, gene) for gene in modified_symbols
            ]

            unmodify_alias_symbols = len(try_modified_symbols) - len(
                set(try_modified_symbols)
            )
            print("以上修改后可能造成重名的基因数: ", unmodify_alias_symbols, ", 不修改")
            adata.var["gene_symbols"] = try_modified_symbols
            # 筛选出存在重复的gene_symbols，且为修改后造成的（与原来名称不同）
            remove_change_genes = list(
                adata.var[
                    adata.var.duplicated(subset=["gene_symbols"], keep=False)
                    & (adata.var["gene_symbols"] != adata.var.index.astype(str))
                ].index
            )
            # 这些基因不修改名称，删除键值对
            if len(remove_change_genes) > 0:
                for key in remove_change_genes:
                    if key in HUGO_alias_symbols_dict:
                        del HUGO_alias_symbols_dict[key]
            # 删除前面添加的 gene_symbols 列
            del adata.var["gene_symbols"]

            modified_symbols = [
                HUGO_alias_symbols_dict.get(gene, gene) for gene in modified_symbols
            ]

            final_duplicate_symbol_count = len(modified_symbols) - len(
                set(modified_symbols)
            )
            print("最终修改后可能造成重名的基因数: ", final_duplicate_symbol_count)
            modify_alias_symbols = modify_alias_symbols - len(remove_change_genes)

        # adata中正式修改基因名
        # ori_index=adata.var.index
        adata.var.index = modified_symbols
        adata.var_names_make_unique()
        # print('最终修改后这些基因发生了变化：',set(ori_index)-set(adata.var.index))
    else:
        modify_alias_symbols = -1
        unmodify_alias_symbols = -1
        final_duplicate_symbol_count = -1
    ## -------------------------------------------------------------------- ##
    # 基因白名单检测，在dict中的物种计算
    if species in check_white_list_species:
        print("——————正在加载基因白名单——————")
        # 读取该物种的基因白名单列表，添加到 gene_white_list
        check_white_list_dir = check_white_list_species[species]
        gene_white_list = []
        for temp_index in range(0, len(check_white_list_dir), 1):
            temp_gene_white_list = pd.read_table(
                check_white_list_dir[temp_index], header=None
            )

            temp_gene_white_list = list(temp_gene_white_list[0])
            gene_white_list.extend(temp_gene_white_list)
            print(
                species,
                "添加第",
                temp_index + 1,
                "个白名单列表后，基因白名单总基因数：",
                len(gene_white_list),
            )

        print("——————进行基因白名单检测——————")
        standard_symbols_count = len(set(modified_symbols) & set(gene_white_list))
        standard_symbols_propotion = round((standard_symbols_count / gene_count), 4)
        print("通过基因白名单检测的基因数占数据集比例：", standard_symbols_propotion)
        unstandard_symbols = ", ".join(set(modified_symbols) - set(gene_white_list))
        print("非白名单的基因名示例：", unstandard_symbols[:1000], "（截取前1000个字符，最后一个基因名不一定完整）")
    else:
        standard_symbols_count = -1
        standard_symbols_propotion = -1

    # 汇总统计
    if "filtered" in list(adata.layers):
        store_raw_count = "True"
    else:
        store_raw_count = "False"

    summary_stat_df.loc[a] = [
        rds_id,
        species,
        store_raw_count,
        cell_count,
        gene_count,
        standard_symbols_count,
        standard_symbols_propotion,
        modify_wrong_format_symbol_count,
        modify_alias_symbols,
        unmodify_alias_symbols,
        final_duplicate_symbol_count,
    ]
    ## -------------------------------------------------------------------- ##
    if 'n_genes' in adata.obs_keys() and 'n_genes_by_counts' in adata.obs_keys():
        del adata.obs['n_genes']

    obs_mapping = {
        "n_genes_by_counts": "nFeature_RNA",  # 基因数
        "n_genes": "nFeature_RNA",  # 基因数
        "total_counts": "nCount_RNA",  # UMI数
        "pct_counts_mito": "percent_mt",  # 线粒体比例
        "Sample ID": "sample",  # 样本名
        "cluster": "raw_cluster",  # 原始分群 cluster（seurat从0开始，默认来源于louvain） raw_cluster（从1开始，来源于louvain或leiden,同时存在louvain列）
    }

    obsm_mapping = {
        "X_harmony": "X_pca_harmony",
    }

    layers_mapping = {  # 不论是 count 还是非 count 的原始矩阵，都保存在raw， 且 adata.raw.X 也有
        "filtered": "raw",
        "abnormal_filtered": "raw",
    }


    adata.obs['cluster'] = adata.obs['cluster'].apply(str)
    
    if adata.X.dtype != "float32":
        adata.X = adata.X.astype(np.float32)

    if "filtered" in list(adata.layers):
        if adata.layers["filtered"].dtype != "float32":
            adata.layers["filtered"] = adata.layers["filtered"].astype(np.float32)
    if "normalised" in list(adata.layers):
        if adata.layers["normalised"].dtype != "float32":
            adata.layers["normalised"] = adata.layers["normalised"].astype(np.float32)

    if "abnormal_filtered" in list(adata.layers):
        if adata.layers["abnormal_filtered"].dtype != "float32":
            adata.layers["abnormal_filtered"] = adata.layers[
                "abnormal_filtered"
            ].astype(np.float32)

    if "abnormal_filtered" in adata.layers:
        adata.layers["raw"] = adata.layers["abnormal_filtered"]
        del adata.layers["abnormal_filtered"]
    elif "filtered" in adata.layers:
        adata.layers["raw"] = adata.layers["filtered"]
        del adata.layers["filtered"]
    else:
        print(
            sys.exit(
                "error: h5ad 中没有存储 adata.layers['abnormal_filtered'] 或 adata.layers['filtered']"
            )
        )

    # h5ad.obs 的列中，需要改的部分名称
    obs_need_change_lst = list(set(adata.obs_keys()) & set(obs_mapping.keys()))
    print(obs_need_change_lst)


    # 修改h5ad.obs
    a_lst = []
    a_value_lst = []  # 记录需要需修改的键对应的值

    for change_obs in adata.obs_keys():
        if change_obs in obs_need_change_lst:
            a_lst.append(change_obs)
            a_value_lst.append(obs_mapping[change_obs])

    # dict(zip(a_lst,a_value_lst))：两个列表生成一个字典
    adata.obs.rename(columns=dict(zip(a_lst, a_value_lst)), inplace=True)
    adata.obs['gname']=adata.obs['sample']
    # 修改obsm
    if "X_harmony" in adata.obsm:
        adata.obsm["X_pca_harmony"] = adata.obsm["X_harmony"]
        del adata.obsm["X_harmony"]
    else:
        pass

    # 删除差异分析结果
    del adata.uns
    ## -------------------------------------------------------------------- ##
    if kingdom == "plants":
        celltype_obs_list = ["annot_sub2", "annot_sub", "annot_full", "annot_article"]
    else:
        celltype_obs_list = ["annot_sub2", "annot_sub", "annot_full"]

    # 最细层级判断
    for ct_level in range(0, len(celltype_obs_list), 1):
        if celltype_obs_list[ct_level] in list(adata.obs.columns.values):
            most_detailed_celltype_level = celltype_obs_list[ct_level]
            print("该数据集最细层级为: ", most_detailed_celltype_level)
            break

    # 将最细层级细胞类型变成 cluster_standard，并转简写
    adata.obs.rename(
        columns=dict(zip([most_detailed_celltype_level], ["cluster_standard"])),
        inplace=True,
    )

    convert_abbr_dict, convert_reverse_dict = convert_abbr()
    print("——————细胞注释全称转简写——————")
    adata.obs["cluster"] = adata.obs["cluster_standard"]
    adata.obs.replace({"cluster": convert_abbr_dict}, inplace=True)

    os.makedirs(outdir, exist_ok=True)
    new_h5ad = outdir + "/" + rds_id + "_lims.h5ad"
    adata.write(filename=new_h5ad, compression="lzf")
    print("——————以下数据集处理完成：", rds_id, " ——————")


summary_stat_df_filename = outdir + "/h5ad_summary_stat.csv"
summary_stat_df.to_csv(summary_stat_df_filename, header=True, index=False)
print("——————输出汇总统计信息，请仔细检查——————")