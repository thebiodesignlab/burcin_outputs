clear all; close all; clc;
mkdir oldAF2

%% Input Sequences

%MBP='KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTRITK'
TMBP='KIEEGKIVFAVGGAPNEIEYWKGVIAEFEKKYPGVTVELKRQATDTEQRRLDLVNALRGKSSDPDVFLMDVAWLGQFIASGWLEPLDDYVQKDNYDLSVFFQSVINLADKQGGKLYALPVYIDAGLLYYRKDLLEKYGYSKPPETWQELVEMAQKIQSGERETNPNFWGFVWQGKQYEGLVCDFVEYVYSNGGSLGEFKDGKWVPTLNKPENVEALQFMVDLIHKYKISPPNTYTEMTEEPVRLMFQQGNAAFERNWPYAWGLHNADDSPVKGKVGVAPLPHFPGHKSAATLGGWHIGISKYSDNKALAWEFVKFVESYSVQKGFAMNLGWNPGRVDVYDDPAVVSKSPHLKELRAVFENAVPRPIVPYYPQLSEIIQKYVNSALAGKISPQEALDKAQKEAEELVKQYS'

cpegfp='SYNVFIMADKQKNGIKANFKIRHNIEDGSVQLAYHYQQNTPIGDGPVLLPDNHYLSVQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKGGTGGSMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYIQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFN'
linker1='VAS';
linker2='AST';

for i=0:length(TMBP)
    MBPN=TMBP(1:i);
    MBPC=TMBP(i+1:end);
    ind=int2str(i);
    foldername=strcat('fusion',ind);
    mkdir(sprintf('oldAF2/%s',foldername));    
    filename=strcat('oldAF2/',foldername,'/fusion',ind,'.fasta');
    fid=fopen(filename,'a');
    fwrite(fid,strcat('>',filename(1:end-6)));
    fwrite(fid,double(sprintf('\n')));
    hybrid{i+1}=strcat(MBPN,linker1,cpegfp,linker2,MBPC);
    fwrite(fid,hybrid{i+1});
    fclose(fid);
    %command=sprintf('alphafold --fasta_paths=%s ----use_gpu=false --max_template_date=2022-01-01 --model_preset=monomer --output_dir=/camp/home/acarb/MBP_EGFP/oldAF2/%s',filename,foldername);
    %status=unix(command)
    clear MBPN MB filename fid foldername
end

quit
