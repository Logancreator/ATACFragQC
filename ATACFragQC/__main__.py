import os, sys, re, getopt, functools, pysam
import pandas as pd
import numpy as np
from plotnine import *
from PIL import Image

from ATACFragQC import __version__

class ArgumentList:
    file_bam = ""
    file_ref = ""
    file_out = False
    delmt = False
    group_marker = ""
    quality = 5
    isize = 147
    cn_len = 10
    chr_filter = ""
    chr_list = ""
    calc_mode = "F"
    pic_list = "a,b,c"
    widthtss = 900
    def __init__(self):
        self.file_bam = ""
        self.file_ref = ""
        self.file_out = False
        self.delmt = False
        self.group_marker = ""
        self.quality = 5
        self.isize = 147
        self.cn_len = 10
        self.calc_mode = "F"
        self.chr_filter = ""
        self.chr_list = ""
        self.pic_list = "a,b,c"
        self.widthtss = 900

def chr_cmp(a, b):
    sa = str(a)
    sb = str(b)
    la = len(sa)
    lb = len(sb)
    lm = min(la, lb)
    for i in range(0, lm):
        if sa[i] != sb[i]:
            oa = ord(sa[i]) if sa[i] != "M" and sa[i] != "m" else 0x7A
            ob = ord(sb[i]) if sb[i] != "M" and sb[i] != "m" else 0x7A
            if oa < 0x3A and oa > 0x2F and ob < 0x3A and ob > 0x2F and la != lb:
                return la - lb
            cd = oa - ob
            return cd
    return la - lb

def bedScan(args):
    (pathname, extension) = os.path.splitext(args.file_bam)
    (filepath, filename) = os.path.split(pathname)
    if not os.path.isfile(args.file_bam+".bai"):
        print("Creating the index file of bam...")
        if os.system("which samtools") == 0:
            os.system("samtools index "+args.file_bam)
        else:
            print("Failed, please check your environment variables...")
            return
    print("Processing the reference...")
    ref_raw = pd.read_table(args.file_ref, comment="#", header=None, dtype={0:str})
    ref_raw.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    type_list = list(set(ref_raw["type"]))
    type_test = "transcript" if "transcript" in type_list else "gene"
    if not type_test in type_list:
        print("There is no suitable term in the gtf to confirm TSSs...")
        return
    chr_list = list(set(ref_raw["seq_id"]))
    if args.chr_list != "":
        chr_list = list(set(chr_list).intersection(set(args.chr_list.split(","))))
    if args.cn_len > 0:
        chr_list = [x for x in chr_list if len(x) < min(args.cn_len, min(list(map(lambda x: len(x), chr_list)))*4)]
    if len(chr_list) == 0:
        print("There is no chromosome would be calculated...")
        return
    fs = pysam.AlignmentFile(args.file_bam, "rb")
    chr_detect = []
    fs_header = fs.header.to_dict()
    if "SQ" in fs_header.keys():
        for term in fs_header["SQ"]:
            if "SN" in term.keys():
                chr_detect.append(str(term["SN"]))
    if len(chr_detect) > 0:
        chr_list = list(set(chr_list).intersection(set(chr_detect)))
    if len(chr_list) == 0:
        print("There is no chromosome would be calculated...")
        return
    chr_list = sorted(list(set(chr_list).difference(set(args.chr_filter.split(",")))), key=functools.cmp_to_key(chr_cmp))
    chr_list_frag = [x for x in chr_list if re.match(r".*[Mm]+,*", x) == None]
    ref = pd.DataFrame()
    for term in chr_list_frag:
        ref = pd.concat([ref, ref_raw[(ref_raw["type"] == type_test) & (ref_raw["strand"] != "-") & ref_raw["seq_id"].str.match("^"+term+"$") & (ref_raw["start"] > args.widthtss+100)]])
    ref = ref[["seq_id", "strand", "start", "end"]].sort_values(by=["seq_id", "start", "strand"])
    ref = ref.drop_duplicates(subset=["seq_id", "start", "strand"], keep="first")
    ref.loc[ref["strand"] == "+", "start"] -= args.widthtss+100
    ref.loc[ref["strand"] == "+", "end"] = ref.loc[ref["strand"] == "+", "start"] + (args.widthtss+100)*2
    ref.index = list(range(ref.index.size))
    
    print("Scaning the distribution of fragments...")
    chr_count = {}
    len_count = [0] * 501
    chr_dist = chr_list_frag if args.delmt else chr_list
    for chr in chr_dist:
        count = 0
        for read in fs.fetch(chr):
            if (read.flag & 65) == 65 and (read.flag & 784) == 0 and read.mapq > args.quality and read.isize > 0 and read.isize < 501:
                len_count[read.isize] += 1
                count += 1
        chr_count[chr] = count
    len_count = pd.DataFrame({"V1": list(range(1, 501)), "V2": len_count[1:]})
    chr_count = pd.DataFrame({"V1": list(chr_count.keys()), "V2": list(chr_count.values())})
    chr_count["V1"] = pd.Categorical(chr_count["V1"], categories=chr_dist, ordered=True)
    
    print("Scaning the fragments around TSSs...")
    dist_count = []
    dist_count_group = {}
    tss_range = (args.widthtss+100)*2+1
    for index, row in ref.iterrows():
        count = np.zeros(tss_range, dtype=np.int64)
        for frag in fs.fetch(row["seq_id"], row["start"], row["end"]):
            if args.calc_mode == "F":
                if (frag.flag & 65) == 65 and (frag.flag & 1808) == 0 and frag.mapq > args.quality and frag.isize > 0 and frag.isize < args.isize:
                    count[range(max(0, frag.pos - row["start"]), min(tss_range, frag.pos + frag.isize - row["start"]))] += 1
            else:
                if (frag.flag & 65) == 65 and (frag.flag & 1808) == 0 and frag.mapq > args.quality and frag.isize > 0:
                    count[max(0, frag.pos - row["start"])] += 1
                    count[min(tss_range, frag.pos + frag.isize - row["start"]) - 1] += 1
        if sum(count) > 0:
            dist_count.append(count)
    dist_count = np.array(dist_count)
    if dist_count.size > 0:
        factors = np.mean(dist_count[:, list(range(0, 100))+list(range(tss_range-100, tss_range))], axis=1)
        factors[factors == 0] = np.mean(factors) if np.mean(factors) > 0 else 1
        dist_count = pd.DataFrame({"V1": list(range(-args.widthtss, args.widthtss+1)), "V2": list(np.mean(dist_count[:, 100:(tss_range-100)] / factors.reshape(len(factors), 1), axis=0))})
    
    print("Saving the results...")
    if args.file_out:
        chr_count.to_csv(pathname+"_chr.tsv", sep="\t", index=False, header=False)
        len_count.to_csv(pathname+"_fl.tsv", sep="\t", index=False, header=False)
        if dist_count.size > 0:
            dist_count.to_csv(pathname+"_tss.tsv", sep="\t", index=False, header=False)
            pd.DataFrame({"V1": factors}).to_csv(pathname+"_base.tsv", sep="\t", index=False, header=False)
    pic_list = ["a", "b", "c"]
    pic_list = list(set(pic_list).intersection(set(args.pic_list.split(","))))
    pic_list.sort()
    if dist_count.size == 0:
        pic_list.remove("c")
    if len(pic_list) == 0:
        return
    if "a" in pic_list:
        ggsave(plot=ggplot(chr_count, aes(x="V1", y="V2"))+geom_bar(stat="identity", width=0.8, fill="#80B1D3")+
            labs(x="Chromosome", y="Fragments")+
            theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), 
            axis_text_y=element_text(colour="black"), axis_text_x=element_text(angle=270, hjust=0.3, vjust=1, colour="black")), 
            width=4, height=6, dpi=200, filename=pathname+".tmpa.png", limitsize=False, verbose=False)
    if "b" in pic_list:
        ggsave(plot=ggplot(len_count, aes(x="V1", y="V2"))+geom_bar(stat="identity", colour="#80B1D3", fill="#80B1D3")+
            labs(x="\nInsert Size", y="Fragments")+
            theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), axis_text=element_text(colour="black")), 
            width=6, height=6, dpi=200, filename=pathname+".tmpb.png", limitsize=False, verbose=False)
    if "c" in pic_list:
        break_length = 50
        if args.widthtss > 6000:
            break_length = 2000
        if args.widthtss > 3000:
            break_length = 1000
        elif args.widthtss > 1000:
            break_length = 500
        elif args.widthtss > 400:
            break_length = 200
        break_range = max(1, args.widthtss//break_length)*break_length
        ggsave(plot=ggplot(dist_count, aes(x="V1", y="V2"))+geom_line(size=1, colour="#80B1D3")+
            labs(x="\nDistance from TSS (bp)", y="Mean TSS enrichment score")+
            scale_x_continuous(breaks=range(-break_range, break_range+1, break_length))+
            theme(plot_title=element_blank(), panel_background=element_blank(), axis_line=element_line(colour="black"), axis_text=element_text(colour="black")), 
            width=6, height=6, dpi=200, filename=pathname+".tmpc.png", limitsize=False, verbose=False)
    width = 0
    height = 0
    imgs = [Image.open(pathname+".tmp"+c+".png") for c in pic_list]
    for img in imgs:
        width += img.width
        height = max(height, img.height)
    result = Image.new("RGBA", (width, height), "white")
    pos = 0
    for img in imgs:
        result.paste(img, (pos, 0))
        pos += img.width
        img.close()
    result.save(pathname+"_qc.png")
    [os.remove(pathname+".tmp"+c+".png") for c in pic_list]

def main():
    opts, args = getopt.getopt(sys.argv[1:], 
        "hodi:r:g:m:q:l:f:c:p:n:w:", 
        ["help", "output", "input=", "reference=", "group=", "mode=", "quality=", "length=", "filter=", "chr=", "pic=", "nl=","widthtss="])
    arguments = ArgumentList()
    help_flag = False
    help_info = "Usage:\nATACFragQC [options] -i <input.bam> -r <reference.gtf>\nArguments:\n"\
        +"-h, --help\t\tShow this help information\n"\
        +"-i, --input <file>\tA aligned & deduped BAM file\n"\
        +"-r, --reference <file>\tGTF genome annotation\n"\
        +"-g, --group [marker]\tThe TSS of each group would be calculated if -g was set\n"\
        +"-o, --output [T/F]\tThe table of results would be saved if -o was set (default: False)\n"\
        +"-m, --mode [F/C]\tThe TSS enrichment would be calculated by fragments (F) or cutsites(C) (default: F)\n"\
        +"-q, --quality [1-255]\tThe quality limit of alignment (default: 5)\n"\
        +"-l, --length [50-500]\tThe length limit of nucleosome-free fragment (default: 147)\n"\
        +"-c, --chr [aaa,bbb]\tThe list of chromosomes would be used (default: all)\n"\
        +"-n, --nl [0-X]\tThe length limit of chromosome names (default: 10, 0 is no limit)\n"\
        +"-p, --pic [a,b,c]\tThe list of images would be shown (default: all)\n"\
        +"-f, --filter [aaa,bbb]\tThe list of chromosomes which should be filtered (default: none)\n"\
        +"-w, --widthtss [100-10000]\tThe width of regions around the TSS (default: 900)\n"\
        +"-d, --delmt\tRemove mtDNA in fragment distribution (default: False)\n"
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help_flag = True
        elif opt in ("-i", "--input"):
            arguments.file_bam = arg
        elif opt in ("-r", "--reference"):
            arguments.file_ref = arg
        elif opt in ("-g", "--group"):
            arguments.group_marker = arg
        elif opt in ("-o", "--output"):
            arguments.file_out = True
        elif opt in ("-d", "--delmt"):
            arguments.delmt = True
        elif opt in ("-q", "--quality"):
            if int(arg) >= 1 and int(arg) <= 255:
                arguments.quality = int(arg)
        elif opt in ("-l", "--length"):
            if int(arg) >= 50 and int(arg) <= 500:
                arguments.isize = int(arg)
        elif opt in ("-m", "--mode"):
            arguments.calc_mode = arg
        elif opt in ("-f", "--filter"):
            arguments.chr_filter = arg
        elif opt in ("-c", "--chr"):
            arguments.chr_filter = arg
        elif opt in ("-p", "--pic"):
            arguments.pic_list = arg
        elif opt in ("-n", "--nl"):
            arguments.cn_len = int(arg)
        elif opt in ("-w", "--widthtss"):
            if int(arg) >= 100 and int(arg) <= 10000:
                arguments.widthtss = int(arg)
    print("ATACFragQC - Version: "+__version__)
    if help_flag or arguments.file_bam == "" or arguments.file_ref == "":
        print(help_info)
        sys.exit()
    bedScan(arguments)

if __name__ == "__main__":
    main()
