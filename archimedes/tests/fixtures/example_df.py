import pandas as pd

def example_df():
    return pd.DataFrame(
        {
            "0": {
                "0":7.3628382683,"1":5.6609773636,"2":6.9171485901,"3":5.2869787216,"4":5.88692379,
                "5":13.3874797821,"6":6.2971229553,"7":3.2141144276,"8":3.2316823006,"9":5.7822804451,
                "10":2.2435512543,"11":1.7055981159,"12":12.8511066437,"13":5.3066930771,"14":11.2680082321,
                "15":2.8846020699,"16":4.6086468697,"17":3.5748836994,"18":5.4393901825,"19":6.0238690376,
                "20":14.2446889877,"21":13.9720697403,"22":3.6030702591,"23":15.5510540009,"24":5.0378332138,
                "25":7.2822446823,"26":6.8934235573,"27":6.5569458008,"28":5.2502856255,"29":14.9126319885,
                "30":3.9893736839,"31":6.6055107117,"32":14.547000885,"33":13.680680275,"34":2.3881351948,
                "35":6.1825251579,"36":3.8028862476,"37":1.1174815893,"38":8.9891490936,"39":15.1695709229,
                "40":5.0741558075,"41":6.3063368797,"42":2.892441988,"43":3.2514247894,"44":3.720095396,
                "45":14.9021911621,"46":4.5512285233,"47":2.1144428253,"48":13.7612485886,"49":2.4658033848,
                "50":4.2058596611,"51":12.932305336,"52":7.5802350044,"53":7.633936882,"54":14.0821418762,
                "55":2.7272927761,"56":15.2532691956,"57":4.5007138252,"58":13.3108205795,"59":6.1096343994,
                "60":6.5757350922,"61":15.6932687759,"62":1.1114404202,"63":16.7076625824,"64":4.2073488235,
                "65":15.4804201126,"66":1.7970142365,"67":14.4915523529,"68":1.932264924,"69":15.1475334167,
                "70":3.2751121521,"71":13.4900989532,"72":14.0973749161,"73":6.5862016678,"74":11.1060638428,
                "75":2.8495998383,"76":15.1968946457,"77":7.4578156471,"78":4.0741829872,"79":4.4401674271,
                "80":14.584438324,"81":2.1347088814,"82":14.012465477,"83":3.1991987228,"84":5.147526741,
                "85":8.679608345,"86":4.8421063423,"87":3.7602100372,"88":4.7240934372,"89":5.93810606,
                "90":2.7895245552,"91":15.3325252533,"92":3.0531189442,"93":16.1812210083,"94":3.8446896076,
                "95":5.2228956223,"96":13.0838499069,"97":6.0543632507,"98":5.7369499207,"99":1.8070656061
                },
            "1": {
                "0":-2.9172041416,"1":10.3659591675,"2":5.5499849319,"3":-1.948782444,"4":-2.1850457191,
                "5":2.2143821716,"6":11.5787858963,"7":7.142967701,"8":11.8473701477,"9":12.7214441299,
                "10":13.0029878616,"11":12.5450468063,"12":4.9667754173,"13":6.651558876,"14":4.1236462593,
                "15":-2.1320414543,"16":7.3765654564,"17":10.8785648346,"18":7.9551939964,"19":-5.3585786819,
                "20":7.4097247124,"21":4.0945067406,"22":12.1024971008,"23":6.8898172379,"24":-4.372071743,
                "25":11.5432109833,"26":12.2258453369,"27":5.7284264565,"28":14.2294158936,"29":1.3325623274,
                "30":4.4748387337,"31":11.3139047623,"32":3.164757967,"33":8.9938907623,"34":-3.2631177902,
                "35":13.4150447845,"36":13.7908859253,"37":11.5255775452,"38":-4.558804512,"39":2.5878703594,
                "40":-0.6143532395,"41":-4.752266407,"42":12.446606636,"43":-3.0963418484,"44":-0.1337213963,
                "45":4.9814143181,"46":0.1962183416,"47":11.800942421,"48":3.1468715668,"49":4.8161253929,
                "50":5.5520987511,"51":3.5104250908,"52":12.4465875626,"53":12.8346204758,"54":5.1203422546,
                "55":13.7832679749,"56":4.7023491859,"57":2.6083626747,"58":3.3607890606,"59":-5.0108604431,
                "60":-3.8125579357,"61":4.9559335709,"62":13.4116258621,"63":5.1544327736,"64":12.5939874649,
                "65":5.1918725967,"66":2.3547198772,"67":2.3923449516,"68":4.9450421333,"69":3.7537620068,
                "70":-4.0196390152,"71":7.3025660515,"72":3.6433463097,"73":12.8949708939,"74":6.3094215393,
                "75":12.3084554672,"76":4.9382696152,"77":14.5019445419,"78":-3.1003260612,"79":11.7104330063,
                "80":6.6284151077,"81":12.288479805,"82":5.2824692726,"83":-1.6829124689,"84":13.3847723007,
                "85":-5.6456141472,"86":12.5845556259,"87":-1.4734385014,"88":13.1881694794,"89":12.0215454102,
                "90":5.934709549,"91":7.2620806694,"92":13.6121492386,"93":5.351995945,"94":11.1434211731,
                "95":11.2054481506,"96":5.3346891403,"97":13.3838176727,"98":-2.5507440567,"99":10.8609962463
                },
            "leiden": {
                "0":"4","1":"0","2":"9","3":"7","4":"4",
                "5":"8","6":"0","7":"5","8":"1","9":"0",
                "10":"1","11":"1","12":"8","13":"9","14":"8",
                "15":"6","16":"10","17":"1","18":"10","19":"4",
                "20":"2","21":"3","22":"1","23":"2","24":"6",
                "25":"0","26":"0","27":"9","28":"0","29":"3",
                "30":"12","31":"0","32":"3","33":"17","34":"6",
                "35":"0","36":"1","37":"1","38":"4","39":"3",
                "40":"7","41":"4","42":"1","43":"6","44":"7",
                "45":"2","46":"7","47":"1","48":"3","49":"5",
                "50":"14","51":"8","52":"0","53":"0","54":"2",
                "55":"1","56":"3","57":"12","58":"8","59":"4",
                "60":"4","61":"2","62":"1","63":"3","64":"0",
                "65":"2","66":"13","67":"3","68":"5","69":"3",
                "70":"6","71":"2","72":"2","73":"0","74":"15",
                "75":"1","76":"2","77":"0","78":"6","79":"0",
                "80":"2","81":"1","82":"2","83":"6","84":"0",
                "85":"4","86":"1","87":"7","88":"0","89":"0",
                "90":"5","91":"2","92":"1","93":"2","94":"1",
                "95":"11","96":"8","97":"0","98":"7","99":"1"
                },
            "sample": {
                "0":"sample4","1":"sample1","2":"sample2","3":"sample2","4":"sample4",
                "5":"sample4","6":"sample5","7":"sample1","8":"sample3","9":"sample1",
                "10":"sample2","11":"sample3","12":"sample1","13":"sample1","14":"sample4",
                "15":"sample3","16":"sample2","17":"sample3","18":"sample5","19":"sample1",
                "20":"sample5","21":"sample4","22":"sample2","23":"sample1","24":"sample5",
                "25":"sample2","26":"sample1","27":"sample1","28":"sample5","29":"sample4",
                "30":"sample5","31":"sample4","32":"sample3","33":"sample5","34":"sample2",
                "35":"sample3","36":"sample5","37":"sample4","38":"sample5","39":"sample3",
                "40":"sample4","41":"sample1","42":"sample2","43":"sample2","44":"sample1",
                "45":"sample2","46":"sample1","47":"sample2","48":"sample4","49":"sample2",
                "50":"sample2","51":"sample2","52":"sample2","53":"sample5","54":"sample4",
                "55":"sample4","56":"sample1","57":"sample4","58":"sample1","59":"sample2",
                "60":"sample5","61":"sample4","62":"sample3","63":"sample4","64":"sample5",
                "65":"sample4","66":"sample2","67":"sample1","68":"sample2","69":"sample2",
                "70":"sample2","71":"sample5","72":"sample5","73":"sample2","74":"sample4",
                "75":"sample2","76":"sample5","77":"sample3","78":"sample2","79":"sample2",
                "80":"sample3","81":"sample2","82":"sample3","83":"sample5","84":"sample2",
                "85":"sample2","86":"sample5","87":"sample3","88":"sample1","89":"sample1",
                "90":"sample1","91":"sample4","92":"sample4","93":"sample3","94":"sample1",
                "95":"sample2","96":"sample5","97":"sample3","98":"sample5","99":"sample5"
                },
            "group": {"0":"group1","1":"group3","2":"group3","3":"group3","4":"group2",
                "5":"group3","6":"group1","7":"group2","8":"group2","9":"group4",
                "10":"group4","11":"group2","12":"group3","13":"group1","14":"group4",
                "15":"group4","16":"group2","17":"group3","18":"group3","19":"group1",
                "20":"group4","21":"group3","22":"group4","23":"group3","24":"group1",
                "25":"group2","26":"group1","27":"group4","28":"group4","29":"group4",
                "30":"group2","31":"group1","32":"group4","33":"group4","34":"group1",
                "35":"group2","36":"group1","37":"group4","38":"group4","39":"group1",
                "40":"group2","41":"group3","42":"group2","43":"group4","44":"group2",
                "45":"group1","46":"group3","47":"group3","48":"group1","49":"group2",
                "50":"group4","51":"group3","52":"group2","53":"group3","54":"group1",
                "55":"group1","56":"group2","57":"group3","58":"group1","59":"group1",
                "60":"group1","61":"group3","62":"group1","63":"group4","64":"group4",
                "65":"group1","66":"group1","67":"group3","68":"group1","69":"group1",
                "70":"group4","71":"group3","72":"group2","73":"group4","74":"group4",
                "75":"group1","76":"group1","77":"group2","78":"group2","79":"group2",
                "80":"group3","81":"group3","82":"group4","83":"group1","84":"group2",
                "85":"group2","86":"group4","87":"group1","88":"group1","89":"group2",
                "90":"group2","91":"group2","92":"group1","93":"group4","94":"group2",
                "95":"group4","96":"group3","97":"group1","98":"group4","99":"group4"
                }
        })