classdef testSolver < matlab.unittest.TestCase
    properties
        TestData
    end
    
    methods (TestMethodSetup)
        function setupOnce(testCase)
            
            
            J = [(-2e-6 -2e-3) 2e-6; ...
                    2e-3      -2e-6];
                
            testCase.TestData.J = J;

            F = [0.2 ;0];
            
            testCase.TestData.F = F;
            
            M0 = [5e4; 5e7];
            
            testCase.TestData.M0 = M0;
            
        end
    end
    
    methods (Test)
        
        function [J, F, M0] = getParams(testCase)
            J = testCase.TestData.J;
            F = testCase.TestData.F;
            M0 = testCase.TestData.M0;
        end
                                    
        function testInverterSymbolic(testCase) 
            
            [J, F, M0] = getParams(testCase);
            tsol = linspace(0, 1e9, 101);
            
            YS = symbolic(tsol, J, F, M0);
            
            odefunc = @(t, y) odefun(t, y, J, F);            
            [ T , Y ] = ode_solver(odefunc, tsol, M0, ...
                                  'tol', 0.1, 'theta', 0.5, ...
                                  'solver', 'myode');
            
            m = max(max(YS' - Y'));
            testCase.assertTrue(m < 6.5e3)
            
        end
                
    end
    
end 


function dy = odefun(t, y, J, F)

y = y(:);
    
dy = J*y + F;

dy = dy(:);

end

function YS = symbolic(tsol, J, F, M0)
    
    recalculate_solution = false;
    
    if recalculate_solution
        syms M1(t) M2(t)
        Y = [M1; M2];
        odes = (diff(Y) == J*Y + F);
        C = (Y(0) == M0);
        [M1(t), M2(t)] = dsolve(odes, C);
        YS = zeros(2, length(tsol));
        for ii=1:length(tsol)
            YS(1, ii) = double(subs(M1(t), t, tsol(ii)));
            YS(2, ii) = double(subs(M2(t), t, tsol(ii)));
        end
    else
        YS_1 = [50000;51036.9761275705;52004.5939672630;52953.0895344805;53882.8407276965;54794.2179772735;55687.5843930490;56563.2959090059;57421.7014250827;58263.1429461821;59087.9557184324;59896.4683627552;60689.0030057947;61465.8754082591;62227.3950907251;62973.8654569570;63705.5839147886;64422.8419946159;65125.9254655483;65815.1144492644;66490.6835316177;67152.9018720370;67802.0333107644;68438.3364739739;69062.0648768138;69673.4670244108;70272.7865108801;70860.2621163772;71436.1279022327;72000.6133042069;72553.9432239008;73096.3381183617;73628.0140879170;74149.1829622729;74660.0523849112;75160.8258958186;75651.7030125804;76132.8793098726;76604.5464973820;77066.8924961874;77520.1015136308;77964.3541167087;78399.8273040138;78826.6945762545;79245.1260053807;79655.2883023436;80057.3448835170;80451.4559358046;80837.7784804624;81216.4664356581;81587.6706777956;81951.5391016267;82308.2166791751;82657.8455174962;83000.5649152953;83336.5114184262;83665.8188742945;83988.6184851843;84305.0388605318;84615.2060681659;84919.2436845359;85217.2728439466;85509.4122868204;85795.7784070061;86076.4852981520;86351.6447991631;86621.3665387600;86885.7579791566;87144.9244588757;87398.9692347172;87647.9935228979;87892.0965393780;88131.3755393904;88365.9258561891;88595.8409390319;88821.2123904120;89042.1300025546;89258.6817931913;89470.9540406286;89679.0313181224;89882.9965275734;90082.9309325575;90278.9141907020;90471.0243854231;90659.3380570355;90843.9302332480;91024.8744590553;91202.2428260404;91376.1060010966;91546.5332545828;91713.5924879225;91877.3502606561;92037.8718169605;92195.2211116426;92349.4608356211;92500.6524409029;92648.8561650673;92794.1310552654;92936.5349917458;93076.1247109151;93212.9558279431];
        YS_2 = [50000000;50988061.9690400;51956646.4608848;52906089.5770662;53836769.5935791;54749057.3108466;55643316.2014538;56519902.5549621;57379165.6198620;58221447.7427198;59047084.5045750;59856404.8546425;60649731.2413715;61427379.7409159;62189660.1830638;62936876.2746808;63669325.7207116;64387300.3427916;65091086.1955141;65780963.6803998;66457207.6576142;67120087.5554768;67769867.4778064;68406806.3091455;69031157.8179041;69643170.7574663;70243088.9652982;70831151.4600973;71407592.5370217;71972641.8610386;72526524.5584262;73069461.3064687;73601668.4213783;74123357.9444794;74634737.7266902;75136011.5113340;75627379.0153148;76109036.0086879;76581174.3926586;77043982.2760393;77497644.0501952;77942340.4625099;78378248.6883979;78805542.4018952;79224391.8448538;79634963.8947697;80037422.1312697;80431926.9012849;80818635.3829355;81197701.6481539;81569276.7240699;81933508.6531822;82290542.5523422;82640520.6705704;82983582.4457315;83319864.5600893;83649500.9947630;83972623.0831080;84289359.5630418;84599836.6283352;84904177.9788905;85202504.8700259;85494936.1607860;85781588.3612971;86062575.6791872;86338010.0650888;86608001.2572416;86872656.8252147;87132082.2127641;87386380.7798434;87635653.8437848;87880000.7196654;88119518.7598762;88354303.3929093;88584448.1613785;88810044.7592882;89031183.0685660;89247951.1948736;89460435.5027096;89668720.6498187;89872889.6209209;90073023.7607745;90269202.8065848;90461504.9197737;90650006.7171197;90834783.3012845;91015908.2907344;91193453.8490720;91367490.7137869;91538088.2244394;91705314.3502864;91869235.7173617;92029917.6350211;92187424.1219625;92341817.9317330;92493160.5777301;92641512.3577106;92786932.3778138;92929478.5761108;93069207.7456881;93206175.5572749];
        YS = (horzcat(YS_1, YS_2))';
    end
    
    
            
end



% M1(t) = exp((t*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025))/23611832414348226068480)*((25*559748704483118570224032497489567596417^(1/2))/1180591620717411303424 - 590295810358705614375/1180591620717411303424)*((262144*(17970445539593863215958248504915676045103853081899175725*559748704483118570224032497489567596417^(1/2) + 425162084916197702796669326240463893893815903049965840768650047017031189523))/(1999102516011137750800116062462741415775*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025)*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025)) - (590295810358705651712*559748704483118570224032497489567596417^(1/2)*exp(-(t*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025))/23611832414348226068480))/(69968588060389821278004062186195949552125*(559748704483118570224032497489567596417^(1/2)/23611832414348226068480 - 4731811215835384205/4722366482869645213696))) + exp(-(t*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025))/23611832414348226068480)*((25*559748704483118570224032497489567596417^(1/2))/1180591620717411303424 + 590295810358705614375/1180591620717411303424)*((128*559748704483118570224032497489567596417^(1/2)*(10434375*559748704483118570224032497489567596417^(1/2) - 10889035741470473694947046763487845298809))/(13993717612077964255600812437239189910425*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025)) + (590295810358705651712*559748704483118570224032497489567596417^(1/2)*exp((t*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025))/23611832414348226068480))/(69968588060389821278004062186195949552125*(559748704483118570224032497489567596417^(1/2)/23611832414348226068480 + 4731811215835384205/4722366482869645213696)));
%         M2(t) = exp((t*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025))/23611832414348226068480)*((262144*(17970445539593863215958248504915676045103853081899175725*559748704483118570224032497489567596417^(1/2) + 425162084916197702796669326240463893893815903049965840768650047017031189523))/(1999102516011137750800116062462741415775*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025)*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025)) - (590295810358705651712*559748704483118570224032497489567596417^(1/2)*exp(-(t*(559748704483118570224032497489567596417^(1/2) - 23659056079176921025))/23611832414348226068480))/(69968588060389821278004062186195949552125*(559748704483118570224032497489567596417^(1/2)/23611832414348226068480 - 4731811215835384205/4722366482869645213696))) - exp(-(t*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025))/23611832414348226068480)*((128*559748704483118570224032497489567596417^(1/2)*(10434375*559748704483118570224032497489567596417^(1/2) - 10889035741470473694947046763487845298809))/(13993717612077964255600812437239189910425*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025)) + (590295810358705651712*559748704483118570224032497489567596417^(1/2)*exp((t*(559748704483118570224032497489567596417^(1/2) + 23659056079176921025))/23611832414348226068480))/(69968588060389821278004062186195949552125*(559748704483118570224032497489567596417^(1/2)/23611832414348226068480 + 4731811215835384205/4722366482869645213696)));

