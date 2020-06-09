
import os


params = [

    #['bp','ets1','ets1-5','C1C1_dist',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','8,13'],
    #['bp','ets1','ets1-5','stretch',7,'GC_bound/GC_free/GA2_free/GT_free','red/lime/pink/blue','-2,2'],
    #['groove','ets1','ets1-5','minorgw',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','0,18'],
    #['groove','ets1','ets1-5','majorgw',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','6,24'],
    #['ABG','ets1','ets1-5','beta',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','0,90'],
    #['ABG','ets1','ets1-5','gamma',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','-180,180'],
    #['ABG','ets1','ets1-5','zeta',7,'GC_bound/GC_free/GA2_free/GA1_free/GT_free/GG1_free/GG2_free','red/lime/pink/purple/blue/orange/yellow','-10,50'],
    #['interaction','ets1','ets1-5','hbonds_core',7,'wt/ga2/ga1/gt/gg1/gg2','lime/pink/purple/blue/orange/yellow','0,16'],
    #['interaction','ets1','ets1-5','hbonds_core_base',7,'wt/ga2/ga1/gt/gg1/gg2','lime/pink/purple/blue/orange/yellow','0,8'],
    #['interaction','ets1','ets1-5','hbonds_core_backbone',7,'wt/ga2/ga1/gt/gg1/gg2','lime/pink/purple/blue/orange/yellow','0,12'],
    #['interaction','ets1','ets1-5','sa_buried',7,'wt/ga2/ga1/gt/gg1/gg2','lime/pink/purple/blue/orange/yellow','0,2000'],

    #['bp','ets1','ets1-4','C1C1_dist',8,'GC_bound/GG1_free/GG2_free/GT_free','red/lime/orange/yellow/blue','8,13'],

    #['groove','ets1','ets1-0','minorgw',6,'GC_bound/GC_free/GT_free/CT_free','red/lime/blue/cyan','0,18'],
    #['groove','ets1','ets1-0','majorgw',6,'GC_bound/GC_free/GT_free/CT_free','red/lime/blue/cyan','6,24'],
    #['ABG','ets1','ets1-0','beta',7,'GC_bound/GC_free/GT_free/CT_free','red/lime/blue/cyan','0,90'],
    #['ABG','ets1','ets1-0','gamma',7,'GC_bound/GC_free/GT_free/CT_free','red/lime/blue/cyan','-180,180'],
    #['ABG','ets1','ets1-0','zeta',7,'GC_bound/GC_free/GT_free/CT_free','red/lime/blue/cyan','-10,50'],
    #['interaction','ets1','ets1-0','hbonds_core',7,'wt/gt/ct','lime/blue/cyan','0,10'],
    #['interaction','ets1','ets1-0','sa_buried',7,'wt/gt/ct','lime/blue/cyan','0,2000'],

    #['groove','ets1','ets1-2','minorgw',7,'AT_bound/AT_free/GC_free/GA2_free/GA1_free','red/lime/green/pink/purple','0,18'],
    #['groove','ets1','ets1-2','majorgw',7,'AT_bound/AT_free/GC_free/GA2_free/GA1_free','red/lime/green/pink/purple','6,24'],
    #['ABG','ets1','ets1-2','beta',7,'AT_bound/AT_free/GC_free/GA2_free/GA1_free','red/lime/green/pink/purple','0,90'],
    #['ABG','ets1','ets1-2','gamma',7,'AT_bound/AT_free/GC_free/GA2_free/GA1_free','red/lime/green/pink/purple','-180,180'],
    #['ABG','ets1','ets1-2','zeta',7,'AT_bound/AT_free/GC_free/GA2_free/GA1_free','red/lime/green/pink/purple','-10,50'],
    #['interaction','ets1','ets1-2','hbonds_core',7,'wt/mut/ga2/ga1','lime/green/pink/purple','0,16'],
    #['interaction','ets1','ets1-2','hbonds_core_base',7,'wt/mut/ga2/ga1','lime/green/pink/purple','0,8'],
    #['interaction','ets1','ets1-2','hbonds_core_backbone',7,'wt/mut/ga2/ga1','lime/green/pink/purple','0,8'],
    #['interaction','ets1','ets1-2','sa_buried',7,'wt/mut/ga2/ga1','lime/green/pink/purple','0,2000'],

    ['bp','p53','p53-1','C1C1_dist',5,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','7,13'],
    ['groove','p53','p53-1','minorgw',3,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','0,18'],
    ['groove','p53','p53-1','majorgw',3,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','6,24'],
    ['ABG','p53','p53-1','beta',5,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','0,90'],
    ['ABG','p53','p53-1','gamma',5,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','-180,180'],
    ['ABG','p53','p53-1','zeta',5,'HG_bound/WC_free/CT_free/TT_free/GT_free','red/lime/cyan/pink/blue','-10,50'],
    #['interaction','p53','p53-1','hbonds_core',5,'wt/ct/tt/gt','lime/cyan/pink/blue','0,40'],
    #['interaction','p53','p53-1','sa_buried',5,'wt/ct/tt/gt','lime/cyan/pink/blue','0,3000'],

    #['groove','myc_max','myc_max','minorgw',8,'AT_bound/AT_free/GT_free/CT_free','red/lime/blue/cyan','0,18'],
    #['groove','myc_max','myc_max','majorgw',8,'AT_bound/AT_free/GT_free/CT_free','red/lime/blue/cyan','6,24'],
    #['ABG','myc_max','myc_max','beta',8,'AT_bound/AT_free/GT_free/CT_free','red/lime/blue/cyan','0,90'],
    #['ABG','myc_max','myc_max','gamma',8,'AT_bound/AT_free/GT_free/CT_free','red/lime/blue/cyan','-180,180'],
    #['ABG','myc_max','myc_max','zeta',8,'AT_bound/AT_free/GT_free/CT_free','red/lime/blue/cyan','-10,50'],
    #['interaction','myc_max','myc_max','hbonds_core',8,'at/gt/ct','lime/blue/cyan','0,30'],
    #['interaction','myc_max','myc_max','sa_buried',8,'at/gt/ct','lime/blue/cyan','1000,2000'],

    #['groove','ctcf','ctcf-1','minorgw',6,'GC_bound/GC_free/GT_free/GG1_free/GG2_free','red/lime/blue/orange/yellow','0,18'],
    #['groove','ctcf','ctcf-1','majorgw',6,'GC_bound/GC_free/GT_free/GG1_free/GG2_free','red/lime/blue/orange/yellow','6,24'],
    #['ABG','ctcf','ctcf-1','beta',6,'GC_bound/GC_free/GT_free/GG1_free/GG2_free','red/lime/blue/orange/yellow','0,90'],
    #['ABG','ctcf','ctcf-1','gamma',6,'GC_bound/GC_free/GT_free/GG1_free/GG2_free','red/lime/blue/orange/yellow','-180,180'],
    #['ABG','ctcf','ctcf-1','zeta',6,'GC_bound/GC_free/GT_free/GG1_free/GG2_free','red/lime/blue/orange/yellow','-10,50'],
    #['interaction','ctcf','ctcf-1','hbonds_core',8,'gc/gt/gg1/gg2','lime/blue/orange/yellow','0,40'],
    #['interaction','ctcf','ctcf-1','sa_buried',8,'gc/gt/gg1/gg2','lime/blue/orange/yellow','1000,3000'],

    #['groove','tbp','tbp','minorgw',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free','red/lime/orange/yellow','0,18'],
    #['groove','tbp','tbp','majorgw',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free','red/lime/orange/yellow','6,24'],
    #['ABG','tbp','tbp','beta',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free','red/lime/orange/yellow','0,90'],
    #['ABG','tbp','tbp','gamma',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free','red/lime/orange/yellow','-180,180'],
    #['ABG','tbp','tbp','zeta',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free/','red/lime/orange/yellow','-10,50'],
    #['interaction','tbp','tbp','hbond_num',9,'AdMLP_bound/AdMLP_free/AC+_free/AC_free','red/lime/orange/yellow','0,5'],

    #['bp','tbp','tbp','C1C1_dist',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','7,13'],
    #['bp','tbp','tbp','C1C1_dist',9,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','7,13'],
    #['groove','tbp','tbp','minorgw',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','0,18'],
    #['groove','tbp','tbp','majorgw',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','6,24'],
    #['ABG','tbp','tbp','beta',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','0,90'],
    #['ABG','tbp','tbp','gamma',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','-180,180'],
    #['ABG','tbp','tbp','zeta',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','-10,50'],
    #['interaction','tbp','tbp','hbond_num',10,'AdMLP_bound/AdMLP_free/CC2_free','red/lime/skyblue','0,5'],
    
]

for i in range(len(params)):
    data_type,tf_name,tf_id,param,resi,mismatch_list,color_list,ylim_list = params[i]
    print("Working on %s %s [%d/%d]"%(tf_id,param,i+1,len(params)))
    if data_type == 'bp':
        os.system("python mkg_MD_violin_bp.py %s %s %s %s %s %s %s"%(tf_name,tf_id,param,resi,mismatch_list,color_list,ylim_list))
    elif data_type == 'groove':
        os.system("python mkg_MD_violin_groove.py %s %s %s %s %s %s %s"%(tf_name,tf_id,param,resi,mismatch_list,color_list,ylim_list))
    elif data_type == 'ABG':
        os.system("python mkg_MD_violin_ABG.py %s %s %s %s %s %s %s"%(tf_name,tf_id,param,resi,mismatch_list,color_list,ylim_list))
    elif data_type == 'interaction':
        os.system("python mkg_MD_violin_interaction.py %s %s %s %s %s %s"%(tf_name,tf_id,param,mismatch_list,color_list,ylim_list))


