clear all
close all
% fname1='1OMP';
% fname2='1ANF';
fname1='4O4B';
%f1=getpdb(fname1);
%%f2=getpdb('4O4B');
%fname2='/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/test_e62db.result/test_e62db_relaxed_rank_5_model_5.pdb';
fname2='/Users/burcinacar/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/test_e62db.result/test_e62db_relaxed_rank_4_model_2.pdb';
%fname2='/Users/acarb/Downloads/S721724_results/model1.pdb';
% [Dist, RMSD, Transf, PDB2TX] = pdbsuperpose(f1, f2,'SEGMENT',{'B', 'A'});
% fname2=PDB2TX;
%pdbwrite('try.pdb',PDB2TX)
chain1='AB';
chain2='A';
[resnum,ndvs,dvs,MSF1f,MSF2f]=difvecPDBseg(fname1,chain1,fname2,chain2);
mode_beg=1;
mode_max=10;
%[MSFi,MSFa,resnum2]=GNM(fname1,[mode_beg:mode_max],chain1);
%ndvs=mine;
for kk=mode_beg:mode_max
    for kkk=mode_beg:mode_max
    ovM(1,kk)=dot(MSF1f(kk,:),ndvs)/(norm(MSF1f(kk,:))*norm(ndvs));
    %ovM(1,kk)=dot(MSF1f(kk,1:370),ndvs)/(norm(MSF1f(kk,1:370))*norm(ndvs));
    ovMsf(kk,kkk)=dot(MSF1f(kk,:),MSF2f(kkk,:))/(norm(MSF1f(kk,:))*norm(MSF2f(kkk,:)));
    end
end

%plot(mode_beg:mode_max,ovM,'LineWidth',6)
bar(mode_beg:mode_max,ovM,'FaceColor',[0.466666666666667 0.674509803921569 0.188235294117647])
set(gca,'FontSize',60)
grid on
axis square
xlabel('Modes of 4O4B')
ylabel('Overlap with Dif Vec')
legend('DMS1')
xlim([mode_beg-0.5 mode_max+0.5])
xtickangle(0)

figure2 = figure('Color',[1 1 1]);
colormap(jet);
axes1 = axes('Parent',figure2);
hold(axes1,'on');
imagesc(ovMsf,[0,1]);
set(gca,'YDir','normal','Fontsize',60);
set(gcf,'Color',[1 1 1])
caxis([0 1]);
axis square
%colormap jet
% Create axes
% Create ylabel
ylabel('4O4B');
% Create xlabel
xlabel('Model');
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 10.5]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.5 10.5]);
box(axes1,'on');
axis(axes1,'square');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'Colormap',...
    [1 1 1;0.999325292432723 0.999290364078894 0.999297005385608;0.998650584865446 0.998580728157788 0.998594010771216;0.997975877298169 0.997871092236681 0.997891016156823;0.997301169730892 0.997161456315575 0.997188021542431;0.996626462163615 0.996451820394469 0.996485026928039;0.995951754596338 0.995742184473363 0.995782032313647;0.995277047029061 0.995032548552257 0.995079037699255;0.994602339461784 0.99432291263115 0.994376043084862;0.993927631894507 0.993613276710044 0.99367304847047;0.993252924327229 0.992903640788938 0.992970053856078;0.992578216759952 0.992194004867832 0.992267059241686;0.991903509192675 0.991484368946726 0.991564064627294;0.991228801625398 0.990774733025619 0.990861070012901;0.990554094058121 0.990065097104513 0.990158075398509;0.989879386490844 0.989355461183407 0.989455080784117;0.989204678923567 0.988645825262301 0.988752086169725;0.98852997135629 0.987936189341195 0.988049091555332;0.987855263789013 0.987226553420088 0.98734609694094;0.987180556221736 0.986516917498982 0.986643102326548;0.986505848654459 0.985807281577876 0.985940107712156;0.985831141087182 0.98509764565677 0.985237113097764;0.985156433519905 0.984388009735664 0.984534118483371;0.984481725952628 0.983678373814557 0.983831123868979;0.983807018385351 0.982968737893451 0.983128129254587;0.983132310818074 0.982259101972345 0.982425134640195;0.982457603250797 0.981549466051239 0.981722140025803;0.98178289568352 0.980839830130133 0.98101914541141;0.981108188116243 0.980130194209027 0.980316150797018;0.980433480548966 0.97942055828792 0.979613156182626;0.979758772981688 0.978710922366814 0.978910161568234;0.979084065414411 0.978001286445708 0.978207166953842;0.978409357847134 0.977291650524602 0.977504172339449;0.977734650279857 0.976582014603496 0.976801177725057;0.97705994271258 0.975872378682389 0.976098183110665;0.976385235145303 0.975162742761283 0.975395188496273;0.975710527578026 0.974453106840177 0.974692193881881;0.975035820010749 0.973743470919071 0.973989199267488;0.974361112443472 0.973033834997965 0.973286204653096;0.973686404876195 0.972324199076858 0.972583210038704;0.973011697308918 0.971614563155752 0.971880215424312;0.972336989741641 0.970904927234646 0.97117722080992;0.971662282174364 0.97019529131354 0.970474226195527;0.970987574607087 0.969485655392434 0.969771231581135;0.97031286703981 0.968776019471327 0.969068236966743;0.969638159472533 0.968066383550221 0.968365242352351;0.968963451905256 0.967356747629115 0.967662247737959;0.968288744337979 0.966647111708009 0.966959253123566;0.967614036770702 0.965937475786903 0.966256258509174;0.966939329203425 0.965227839865796 0.965553263894782;0.966264621636147 0.96451820394469 0.96485026928039;0.96558991406887 0.963808568023584 0.964147274665998;0.964915206501593 0.963098932102478 0.963444280051605;0.964240498934316 0.962389296181372 0.962741285437213;0.963565791367039 0.961679660260265 0.962038290822821;0.962891083799762 0.960970024339159 0.961335296208429;0.962216376232485 0.960260388418053 0.960632301594037;0.961541668665208 0.959550752496947 0.959929306979644;0.960866961097931 0.958841116575841 0.959226312365252;0.960192253530654 0.958131480654735 0.95852331775086;0.959517545963377 0.957421844733628 0.957820323136468;0.9588428383961 0.956712208812522 0.957117328522076;0.958168130828823 0.956002572891416 0.956414333907683;0.957493423261546 0.95529293697031 0.955711339293291;0.956818715694269 0.954583301049204 0.955008344678899;0.956144008126992 0.953873665128097 0.954305350064507;0.955469300559715 0.953164029206991 0.953602355450115;0.954794592992438 0.952454393285885 0.952899360835722;0.954119885425161 0.951744757364779 0.95219636622133;0.953445177857884 0.951035121443673 0.951493371606938;0.952770470290606 0.950325485522566 0.950790376992546;0.952095762723329 0.94961584960146 0.950087382378153;0.951421055156052 0.948906213680354 0.949384387763761;0.950746347588775 0.948196577759248 0.948681393149369;0.950071640021498 0.947486941838142 0.947978398534977;0.949396932454221 0.946777305917035 0.947275403920585;0.948722224886944 0.946067669995929 0.946572409306192;0.948047517319667 0.945358034074823 0.9458694146918;0.94737280975239 0.944648398153717 0.945166420077408;0.946698102185113 0.943938762232611 0.944463425463016;0.946023394617836 0.943229126311504 0.943760430848624;0.945348687050559 0.942519490390398 0.943057436234231;0.944673979483282 0.941809854469292 0.942354441619839;0.943999271916005 0.941100218548186 0.941651447005447;0.943324564348728 0.94039058262708 0.940948452391055;0.942649856781451 0.939680946705973 0.940245457776663;0.941975149214174 0.938971310784867 0.93954246316227;0.941300441646897 0.938261674863761 0.938839468547878;0.94062573407962 0.937552038942655 0.938136473933486;0.939951026512343 0.936842403021549 0.937433479319094;0.939276318945065 0.936132767100442 0.936730484704702;0.938601611377788 0.935423131179336 0.936027490090309;0.937926903810511 0.93471349525823 0.935324495475917;0.937252196243234 0.934003859337124 0.934621500861525;0.936577488675957 0.933294223416018 0.933918506247133;0.93590278110868 0.932584587494911 0.933215511632741;0.935228073541403 0.931874951573805 0.932512517018348;0.934553365974126 0.931165315652699 0.931809522403956;0.933878658406849 0.930455679731593 0.931106527789564;0.933203950839572 0.929746043810487 0.930403533175172;0.932529243272295 0.92903640788938 0.929700538560779;0.931854535705018 0.928326771968274 0.928997543946387;0.931179828137741 0.927617136047168 0.928294549331995;0.930505120570464 0.926907500126062 0.927591554717603;0.929830413003187 0.926197864204956 0.926888560103211;0.92915570543591 0.925488228283849 0.926185565488818;0.928480997868633 0.924778592362743 0.925482570874426;0.927806290301356 0.924068956441637 0.924779576260034;0.927131582734078 0.923359320520531 0.924076581645642;0.926456875166801 0.922649684599425 0.92337358703125;0.925782167599524 0.921940048678318 0.922670592416857;0.925107460032247 0.921230412757212 0.921967597802465;0.92443275246497 0.920520776836106 0.921264603188073;0.923758044897693 0.919811140915 0.920561608573681;0.923083337330416 0.919101504993894 0.919858613959289;0.922408629763139 0.918391869072788 0.919155619344896;0.921733922195862 0.917682233151681 0.918452624730504;0.921059214628585 0.916972597230575 0.917749630116112;0.920384507061308 0.916262961309469 0.91704663550172;0.919709799494031 0.915553325388363 0.916343640887328;0.919035091926754 0.914843689467257 0.915640646272935;0.918360384359477 0.91413405354615 0.914937651658543;0.9176856767922 0.913424417625044 0.914234657044151;0.917010969224923 0.912714781703938 0.913531662429759;0.916336261657646 0.912005145782832 0.912828667815367;0.915661554090369 0.911295509861726 0.912125673200974;0.914986846523092 0.910585873940619 0.911422678586582;0.91281868583413 0.904135063852312 0.905786175074207;0.910650525145169 0.897684253764005 0.900149671561832;0.908482364456207 0.891233443675698 0.894513168049457;0.906314203767246 0.884782633587392 0.888876664537082;0.904146043078284 0.878331823499085 0.883240161024707;0.901977882389323 0.871881013410778 0.877603657512332;0.899809721700361 0.865430203322471 0.871967153999957;0.8976415610114 0.858979393234164 0.866330650487582;0.895473400322438 0.852528583145857 0.860694146975206;0.893305239633477 0.84607777305755 0.855057643462831;0.891137078944515 0.839626962969243 0.849421139950456;0.888968918255554 0.833176152880936 0.843784636438081;0.886800757566592 0.826725342792629 0.838148132925706;0.884632596877631 0.820274532704322 0.832511629413331;0.882464436188669 0.813823722616015 0.826875125900956;0.880296275499708 0.807372912527708 0.821238622388581;0.878128114810746 0.800922102439401 0.815602118876206;0.875959954121785 0.794471292351094 0.809965615363831;0.873791793432823 0.788020482262787 0.804329111851456;0.871623632743862 0.78156967217448 0.798692608339081;0.8694554720549 0.775118862086173 0.793056104826705;0.867287311365939 0.768668051997866 0.78741960131433;0.865119150676977 0.762217241909559 0.781783097801955;0.862950989988016 0.755766431821252 0.77614659428958;0.860782829299054 0.749315621732945 0.770510090777205;0.858614668610093 0.742864811644638 0.76487358726483;0.856446507921131 0.736414001556331 0.759237083752455;0.85427834723217 0.729963191468024 0.75360058024008;0.852110186543208 0.723512381379717 0.747964076727705;0.849942025854247 0.71706157129141 0.74232757321533;0.847773865165285 0.710610761203103 0.736691069702955;0.845605704476324 0.704159951114796 0.731054566190579;0.843437543787362 0.697709141026489 0.725418062678204;0.841269383098401 0.691258330938182 0.719781559165829;0.839101222409439 0.684807520849875 0.714145055653454;0.836933061720478 0.678356710761568 0.708508552141079;0.834764901031516 0.671905900673261 0.702872048628704;0.832596740342555 0.665455090584954 0.697235545116329;0.830428579653593 0.659004280496647 0.691599041603954;0.828260418964632 0.65255347040834 0.685962538091579;0.82609225827567 0.646102660320033 0.680326034579204;0.823924097586709 0.639651850231726 0.674689531066829;0.821755936897747 0.633201040143419 0.669053027554454;0.819587776208786 0.626750230055113 0.663416524042078;0.817419615519824 0.620299419966805 0.657780020529703;0.815251454830863 0.613848609878499 0.652143517017328;0.813083294141901 0.607397799790192 0.646507013504953;0.81091513345294 0.600946989701885 0.640870509992578;0.808746972763978 0.594496179613578 0.635234006480203;0.806578812075017 0.588045369525271 0.629597502967828;0.804410651386055 0.581594559436964 0.623960999455453;0.802242490697094 0.575143749348657 0.618324495943078;0.800074330008132 0.56869293926035 0.612687992430703;0.797906169319171 0.562242129172043 0.607051488918328;0.795738008630209 0.555791319083736 0.601414985405953;0.793569847941248 0.549340508995429 0.595778481893577;0.791401687252286 0.542889698907122 0.590141978381202;0.789233526563325 0.536438888818815 0.584505474868827;0.787065365874363 0.529988078730508 0.578868971356452;0.784897205185402 0.523537268642201 0.573232467844077;0.78272904449644 0.517086458553894 0.567595964331702;0.780560883807479 0.510635648465587 0.561959460819327;0.778392723118517 0.50418483837728 0.556322957306952;0.776224562429556 0.497734028288973 0.550686453794577;0.774056401740594 0.491283218200666 0.545049950282202;0.771888241051633 0.484832408112359 0.539413446769826;0.769720080362671 0.478381598024052 0.533776943257451;0.76755191967371 0.471930787935745 0.528140439745076;0.765383758984749 0.465479977847438 0.522503936232701;0.763215598295787 0.459029167759131 0.516867432720326;0.761047437606825 0.452578357670824 0.511230929207951;0.758879276917864 0.446127547582517 0.505594425695576;0.756711116228902 0.43967673749421 0.499957922183201;0.754542955539941 0.433225927405903 0.494321418670826;0.752374794850979 0.426775117317596 0.488684915158451;0.750206634162018 0.420324307229289 0.483048411646076;0.748038473473057 0.413873497140982 0.477411908133701;0.745870312784095 0.407422687052675 0.471775404621326;0.743702152095133 0.400971876964368 0.46613890110895;0.741533991406172 0.394521066876061 0.460502397596575;0.739365830717211 0.388070256787754 0.4548658940842;0.737197670028249 0.381619446699447 0.449229390571825;0.735029509339288 0.375168636611141 0.44359288705945;0.732861348650326 0.368717826522833 0.437956383547075;0.730693187961365 0.362267016434527 0.4323198800347;0.728525027272403 0.35581620634622 0.426683376522325;0.726356866583442 0.349365396257913 0.42104687300995;0.72418870589448 0.342914586169606 0.415410369497575;0.722020545205519 0.336463776081299 0.409773865985199;0.719852384516557 0.330012965992992 0.404137362472824;0.717684223827596 0.323562155904685 0.398500858960449;0.715516063138634 0.317111345816378 0.392864355448074;0.713347902449673 0.310660535728071 0.387227851935699;0.711179741760711 0.304209725639764 0.381591348423324;0.70901158107175 0.297758915551457 0.375954844910949;0.706843420382788 0.29130810546315 0.370318341398574;0.704675259693827 0.284857295374843 0.364681837886199;0.702507099004865 0.278406485286536 0.359045334373824;0.700338938315904 0.271955675198229 0.353408830861449;0.698170777626942 0.265504865109922 0.347772327349074;0.696002616937981 0.259054055021615 0.342135823836698;0.693834456249019 0.252603244933308 0.336499320324323;0.691666295560058 0.246152434845001 0.330862816811948;0.689498134871096 0.239701624756694 0.325226313299573;0.687329974182135 0.233250814668387 0.319589809787198;0.685161813493173 0.22680000458008 0.313953306274823;0.682993652804212 0.220349194491773 0.308316802762448;0.68082549211525 0.213898384403466 0.302680299250073;0.678657331426289 0.207447574315159 0.297043795737698;0.676489170737327 0.200996764226852 0.291407292225323;0.674321010048366 0.194545954138545 0.285770788712948;0.672152849359404 0.188095144050238 0.280134285200573;0.669984688670443 0.181644333961931 0.274497781688197;0.667816527981481 0.175193523873624 0.268861278175822;0.66564836729252 0.168742713785317 0.263224774663447;0.663480206603558 0.16229190369701 0.257588271151072;0.661312045914597 0.155841093608703 0.251951767638697;0.659143885225635 0.149390283520396 0.246315264126322;0.656975724536674 0.142939473432089 0.240678760613947;0.654807563847712 0.136488663343782 0.235042257101572;0.652639403158751 0.130037853255475 0.229405753589197;0.650471242469789 0.123587043167169 0.223769250076822;0.648303081780828 0.117136233078861 0.218132746564447;0.646134921091866 0.110685422990555 0.212496243052071;0.643966760402905 0.104234612902248 0.206859739539696;0.641798599713943 0.0977838028139406 0.201223236027321;0.639630439024982 0.0913329927256337 0.195586732514946;0.63746227833602 0.0848821826373266 0.189950229002571;0.635294117647059 0.0784313725490197 0.184313725490196],...
    'FontSize',60,'Layer','top');
% Create colorbar
colorbar(axes1);
xlim([mode_beg-0.5 mode_max+.5])
ylim([mode_beg-0.5 mode_max+.5])



