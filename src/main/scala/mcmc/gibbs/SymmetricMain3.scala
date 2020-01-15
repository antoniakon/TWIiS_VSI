package mcmc.gibbs

import breeze.linalg.{*, DenseMatrix, DenseVector, max}
import breeze.numerics.{exp, log, pow, sqrt}
import breeze.stats.mean
import structure.DVStructure

class SymmetricMain3 extends VariableSelection {
  val dv = Array(2.6017059164528966, -0.009771336315506675, 0.007912232574264861, 0.04007021825945167, 2.7441230878433105, 0.027389556986746822, -0.06697554337599647, -8.007114122174142E-4, 0.0, 0.026900437913421865, 3.6658784094663885, 2.49522684662327, -0.018418121343520145, 0.0, 0.02518103384371097, 2.23961135155839, -1.3927526377103536, 0.013779417856511118, 0.03744812080378675, 0.008686020641730074, -0.006342141327408062, 1.783061952079203, -0.015875738391452312, -0.011220586349821161, 0.016730791239580316, 0.0031714521057734893, 0.023208982199357388, 5.578001433242671E-4, 0.0, 3.5629007719930456E-4, 0.8478035636804637, 0.0035385835872639004, 0.03448048158635235, 0.0, 0.016910648617650063, -0.009368709234820406, 0.0030404012698019172, 0.03203897638678914, 0.02562109914990373, 0.006085918534248676, -0.6127898353284323, 0.012286933096967727, -0.07782555618127059, 0.05249466874041374, -0.010116152438496472, -0.014339608855101707, -0.12918199343877793, 0.011483879096421803, 0.0, 0.01288027780692553, 0.0038239283559973352, 0.8405081720024301, -0.028329312840997076, 0.0, -0.01657004966660067, 0.5271052866505676, -0.02400369455385493, -1.5414120944729437, 0.22434385344355495, -0.048651255365600846, 0.03510137768962216, 1.3956769394733277, -0.006208417984605372, 0.004571603754511012, 0.01920259618589719, -0.19731872423386332, 3.0951031257523773, 1.6911727990400254, 0.0, 1.9134228442030587, 0.013617223631881427, -0.005550190903711253, 0.9385665761220506, 0.0, 1.222982841811133, 0.014999044811879131, -0.009671830288291224, -0.021029616651297656, 0.007433814039348221, 0.0022368292335184013, -0.013136232312181284, 0.017199121035394296, -0.01740071491815674, 0.12037636268847685, 0.026776989952518104, 0.011875981660290195, 0.05464804932128788, 0.015713945271984244, 0.0, 0.003320441421002656, 0.021911422466208946, 0.03214958980778366, 0.019004472114151256, 0.0, -0.010538127567225223, -0.01695009252238207, -0.023956683183894063, -0.04021795491570351, 4.117778235193302, 6.638783033116403E-4, 1.6955349206383794, -0.004071167833090188, 0.0013643548981905154, 0.006951240711676783, 0.013893302554700424, 0.002864385926936876, 0.0047829099571852175, 0.020544066082143313, 0.0, 1.2633265516121868, -0.029402884878882516, 0.017648366212642333, -0.015391970610614968, 0.0, 1.3488416793342801, 0.050658742453676166, -0.016220645016616363, 1.3391089108306689, -0.007742778445585529, -0.6136648496914848, 0.0056994343438359595, 0.029390371361697067, -0.10059827247469787, 0.08447235068592644, 8.495069933195136E-4, -0.01933678877529636, 0.03507937954581327, -0.0034912196549974506, 0.0, 0.007013368310608085, 0.08774515792232494, -0.03675755669828546, 0.8197278883875869, 0.0, -0.006347559265661571, 0.0024170448104629885, 0.004621910593119889, 0.004212811224493226, 0.03299560368660157, -0.02051391035057748, -2.1935780084758583, -0.04535858073656455, -0.006759353253553382, -2.69981923138738, -0.01081541402145468, 0.021084648461874776, 1.563805470155378, 3.82602974254243, 0.0, -0.009409404225677142, 0.005196988178242716, 1.568745530604427, 0.01696555599522003, 0.0, -0.01685537612434539, 0.005766787845331351, 1.7747857586180344, -0.03416732272136113, -0.003704758060042889, 0.012454476706899842, -0.028265555962092483, 0.6034788005595431, 0.008340428642646398, -0.039393489793803756, 0.032953553906615835, 0.0138137701765469, -0.010068371380224551, 3.1628248096308025, 0.0, 0.010593599464334882, 0.09156384064108444, -0.041834120824594714, -0.019919418633226764, 0.0, 0.007454394397134417, -0.02778478917093173, -1.4028707000422551, -0.0021023720668093414, -0.01883884516040393, 0.011143166051209684, 0.06675041553683976, 0.04183130383762347, 3.367854295207169E-5, -0.006206358257623055, -0.0023250894466018864, 4.296511236637127, 3.901439261002965, 1.850037666924977, 0.0, -5.114498313832163E-4, 1.070802633175997, -0.053443193597903046, 0.011830425256937291, 0.0, 0.033005551868997715, -0.02714971195045028, -0.021166810819407522, -0.033705053568302536, 5.267473237559164, 1.8518981921016444, 2.5724162689874754, 0.0077822099332529185, -0.9530861427055595, 1.6873965052815756, -0.03926149566915929, 0.007772020547010697, -0.12430173131115425, 7.084771275308971E-4, 0.0, 0.003925500419002989, -0.011227184979535794, -0.05167512036554012, -0.026674460286618445, 0.0, 0.019475567060499282, -0.7179743021041293, -0.004344596906649192, 0.007249345579431036, 0.02133883570545327, -0.0045051312354994575, 0.0031850346480215394, -0.016008757447594613, 0.00659878165759745, -0.01786215656375753, -0.025601715868604387, -0.025482591750035144, 0.01660406697533271, 2.7534886555703064, 0.0, -0.010400119057408695, 9.973076548113538E-4, -0.011153367910363806, -0.039370026798710565, 0.0, -0.008060907741521359, 0.027515018504937698, -0.029355038784475165, 0.030778345612772474, -0.008952523183035173, 0.005633259474114694, 0.011210897476641387, -0.016973193828707055, 2.430669517708118, 0.08191611268617145, 3.0374209166848276, 0.0421678517663094, -0.009039518606554359, 2.4954346434120573, 0.0, 0.018230749251786707, -0.005889079856932727, -0.014872292222471458, 0.17006026830098006, 0.0, 0.009841531368490268, -0.029893817460370148, 0.01355068391733152, 0.017449719203466524, -1.003980247598381, -0.24803454731726943, -0.02840920972623035, 1.2757721588734565, 0.0423466833480538, -0.017744672059155885, 0.024415063181221663, 0.1549691601501992, 0.3937792775610309, 0.006337340324436993, 0.0, -0.07527158715313358, -0.00817329287120582, -0.012788790025210224, 0.003921589933257613, 0.0, -0.03962543698819885, 0.08460029491203659, -0.033204319757496596, 0.027051195646422152, 0.0024014977190241547, 0.07384643276185794, -0.012293786729088539, 0.10703591747743424, 1.9086761465634858, 5.448874234780498, 1.0027838965759301, -0.02158536583495425, -0.017528397924997767, -0.009182908948964041, 0.0, 2.230572420345381, -0.0013334511999952994, 0.03544239860477619, 3.726115133285349, 0.0, 0.023608054013079114, -0.007140256725257687, 2.3146189561417323, -0.6726027254442168, 3.769250237928747, -0.03956786646998975, 2.7187635902268363, 0.800566755732455, 0.019341082575404404, -0.11967774356450467, 0.009787691930742031, -0.004268467741600024, 0.06121970156937635, -0.0020537136987004213, 0.0, -0.006819586103552202, 0.010302669340880997, 0.055360405426384665, -0.054956037788458104, 0.0, 0.0021467449388198724, -0.013874707347192636, 0.021860488383471613, 0.00899292972194786, 0.025003861478492217, 0.023789621029437276, -0.025169323262606378, -0.002598168816693438, -0.01076183623897745, -0.006279324726378747, -0.033855531998819806, 0.010343179590247865, -0.02860605920974769, -0.1780880719722548, 0.0, 0.1255623210685703, -0.024419660658919443, -8.506460577821602E-4, 0.045249311905379876, 0.0, 0.03148124589267364, -0.01800091710266011, 0.057256595358706724, -0.038618132337551646, 0.0202446150655338, 0.03453189965027892, -0.025123214990926433, -0.013609797777575088, -0.016951558099182914, 1.1285450718947145, -0.010310793550872352, 1.476054308391575, 0.660087084973631, -0.030727376136287114, 0.0, -0.018029026140666753, -0.013140143732990455, -0.018853259025706243, -0.01494750360911966, 0.0, 0.6377027312919384, -0.02825611622857848, 0.007640191579344218, 0.03048020620412202, -0.03884596771176862, 0.011821571544450835, 0.013744397445671993, 0.04269355903476545, -0.14765063133071757, -0.019796381142077452, 0.02103176508153761, -0.02168481380165882, -0.0013885111718701481, 0.03944751359931696, 0.0, 0.021131136367325956, 1.139627079186505, -0.013120667737221934, -0.019206570561545143, 0.0, 0.0021449168340867647, 0.012961133070772214, 0.0018698963250457924, 0.008072007383928671, 0.017198628299538474, 0.02901690441255028, 1.2844650324816175E-4, -0.015748168609363053, 0.03285079421632444, 0.008392400614867538, -1.9861004663646382, 1.8149364603160316, -0.012694529266456409, -0.010822860402497858, 0.0, 0.03825187725079332, 2.243292548429352, -0.01851044828307917, 0.050796718498430864, 0.0, 4.337090054433234, -0.0691467126137145, 0.02369826743589345, -1.6747681193944786, 0.09328315951341572, 0.0088103895589982)
  val dm = new DenseMatrix(20, 20, dv)

  val ind = Array(0.9962003799620035, 0.02339766023397669, 0.039396060393960666, 0.04229577042295774, 0.9963003699630034, 0.03469653034696527, 0.1268873112688735, 0.02779722027797215, 0.0, 0.02569743025697434, 0.9970002999700031, 0.9963003699630035, 0.04689531046895304, 0.0, 0.033296670332966787, 0.9966003399660035, 0.9969003099690033, 0.02159784021597842, 0.036996300369963105, 0.02119788021197883, 0.02079792020797921, 0.996100389961004, 0.07269273072692732, 0.023897610238976144, 0.025697430256974223, 0.03409659034096584, 0.025897410258974074, 0.021097890210978864, 0.0, 0.02969703029697031, 0.9905009499050096, 0.057094290570943056, 0.025097490250974897, 0.0, 0.026797320267973185, 0.03789621037896212, 0.05369463053694628, 0.025697430256974237, 0.02489751024897506, 0.027097290270972765, 0.9305069493050674, 0.03359664033596638, 0.15448455154484578, 0.08939106089391063, 0.025497450254974418, 0.02809719028097189, 0.2947705229477048, 0.02869713028697139, 0.0, 0.023397660233976676, 0.06979302069793024, 0.9945005499450094, 0.031196880311968742, 0.0, 0.02969703029697022, 0.7553244675532451, 0.03819618038196164, 0.9970002999700026, 0.410958904109589, 0.0558944105589441, 0.028897110288971014, 0.996800319968003, 0.030296970302969763, 0.025797420257974237, 0.027497250274972566, 0.41005899410058905, 0.9965003499650035, 0.996900309969003, 0.0, 0.9962003799620036, 0.07679232076792321, 0.024097590240975877, 0.9961003899610038, 0.0, 0.9964003599640038, 0.03319668033196681, 0.026497350264973588, 0.027197280271972865, 0.03209679032096772, 0.04169583041695843, 0.021297870212978656, 0.025197480251974838, 0.05169483051694834, 0.2882711728827119, 0.03039696030396967, 0.022997700229977016, 0.024797520247975213, 0.027697230276972365, 0.0, 0.022297770222977666, 0.0682931706829318, 0.02549745025497456, 0.031296870312968746, 0.0, 0.0441955804419559, 0.04119588041195899, 0.025697430256974373, 0.024697530246975224, 0.9966003399660035, 0.0341965803419658, 0.9963003699630041, 0.059194080591940944, 0.02269773022697733, 0.026297370262973754, 0.024897510248975133, 0.0760923907609236, 0.03189681031896813, 0.03409659034096579, 0.0, 0.9961003899610038, 0.021997800219977954, 0.02149785021497852, 0.02329767023297674, 0.0, 0.9963003699630038, 0.03279672032796739, 0.058394160583941604, 0.9965003499650035, 0.02739726027397258, 0.8125187481251887, 0.0416958304169582, 0.024797520247975168, 0.1602839716028397, 0.19688031196880287, 0.022297770222977693, 0.02729727027297274, 0.0317968203179684, 0.044895510448955255, 0.0, 0.025597440255974407, 0.17798220177982202, 0.019398060193980552, 0.9940005999400088, 0.0, 0.02489751024897513, 0.02399760023997591, 0.02509749025097484, 0.01969803019698031, 0.03539646035396466, 0.022597740225977412, 0.9958004199580042, 0.10398960103989581, 0.02089791020897915, 0.9964003599640041, 0.03689631036896304, 0.05489451054894528, 0.996500349965004, 0.9966003399660035, 0.0, 0.02249775022497757, 0.04129587041295845, 0.9956004399560043, 0.028097190280971795, 0.0, 0.027797220277972254, 0.07859214078592151, 0.9965003499650036, 0.03309669033096684, 0.03099690030996904, 0.09769023097690224, 0.025997400259973966, 0.8704129587041306, 0.04009599040096005, 0.07369263073692654, 0.0678932106789324, 0.07859214078592201, 0.043495650434956436, 0.9963003699630036, 0.0, 0.032796720327967065, 0.1867813218678134, 0.0503949605039497, 0.0933906609339068, 0.0, 0.0210978902109789, 0.02609739026097381, 0.9963003699630038, 0.031496850314968496, 0.060193980601939874, 0.02549745025497464, 0.13928607139286067, 0.02429757024297577, 0.07599240075992382, 0.025297470252974744, 0.030496950304969593, 0.9968003199680032, 0.9962003799620042, 0.996200379962004, 0.0, 0.02979702029797025, 0.9962003799620035, 0.023897610238976092, 0.02669733026697335, 0.0, 0.091890810918908, 0.029297070292970764, 0.035496450354964425, 0.03209679032096798, 0.9970002999700028, 0.9969003099690029, 0.9962003799620036, 0.04499550044995495, 0.9958004199580057, 0.9965003499650035, 0.08279172082791736, 0.0337966203379662, 0.2795720427957199, 0.019698030196980285, 0.0, 0.02719728027197285, 0.028497150284971507, 0.09039096090390955, 0.028397160283971597, 0.0, 0.022597740225977263, 0.9859014098590149, 0.024397560243975464, 0.035796420357964084, 0.026697330266973334, 0.028697130286971344, 0.05289471052894689, 0.024297570242975793, 0.025697430256974366, 0.02179782021797824, 0.024997500249975022, 0.03159684031596833, 0.02239776022397758, 0.9964003599640034, 0.0, 0.06959304069593032, 0.018898110188981056, 0.02699730026997309, 0.03499650034996489, 0.0, 0.06249375062493759, 0.02669733026697325, 0.06509349065093464, 0.052194780521947745, 0.03149685031496859, 0.03369663033696625, 0.024097590240975863, 0.023197680231976856, 0.9966003399660035, 0.14398560143985586, 0.9967003299670036, 0.023397660233976575, 0.03259674032596735, 0.9963003699630035, 0.0, 0.0670932906709329, 0.0234976502349765, 0.024697530246975286, 0.3164683531646825, 0.0, 0.02529747025297481, 0.034696530346965, 0.038796120387961104, 0.028097190280971993, 0.9966003399660033, 0.45475452454754595, 0.06589341065893406, 0.9965003499650038, 0.025997400259974046, 0.03259674032596714, 0.046295370462953504, 0.30186981301869803, 0.6205379462053817, 0.0395960403959605, 0.0, 0.1399860013998606, 0.028497150284971465, 0.04319568043195679, 0.0238976102389761, 0.0, 0.0402959704029599, 0.18118188181181957, 0.02999700029996997, 0.025097490250974908, 0.029097090290970917, 0.24297570242975708, 0.038096190380961986, 0.24657534246575322, 0.996100389961004, 0.9969003099690034, 0.9914008599140088, 0.0304969503049695, 0.02559744025597433, 0.03159684031596836, 0.0, 0.9966003399660031, 0.021997800219978006, 0.0395960403959605, 0.996300369963004, 0.0, 0.02419758024197581, 0.021797820217978194, 0.9964003599640037, 0.9645035496450346, 0.9965003499650036, 0.03699630036996309, 0.9963003699630038, 0.9884011598840085, 0.03229677032296771, 0.1804819518048196, 0.020297970202979795, 0.02279772022797723, 0.1399860013998604, 0.034096590340965925, 0.0, 0.02259774022597746, 0.039796020397960104, 0.18638136186381368, 0.02619738026197377, 0.0, 0.0271972802719728, 0.027697230276972237, 0.029397060293970618, 0.04239576042395764, 0.07259274072592749, 0.024297570242975755, 0.025197480251974796, 0.026197380261973754, 0.02559744025597446, 0.02909709029097085, 0.02659734026597344, 0.03499650034996502, 0.033596640335966425, 0.36796320367963214, 0.0, 0.23927607239276147, 0.026597340265973317, 0.024997500249975046, 0.059894010598940214, 0.0, 0.029097090290970937, 0.04659534046595328, 0.0940905909409057, 0.02419758024197568, 0.026797320267973167, 0.037996200379962014, 0.02539746025397451, 0.024597540245975397, 0.042195780421957645, 0.996300369963004, 0.027197280271972837, 0.9960003999600039, 0.9169083091690828, 0.02999700029996993, 0.0, 0.04139586041395863, 0.04689531046895288, 0.032196780321967794, 0.029097090290970948, 0.0, 0.8150184981501838, 0.1479852014798519, 0.02789721027897206, 0.021397860213978624, 0.04539546045395442, 0.0270972902709729, 0.023897610238975967, 0.024597540245975387, 0.3212678732126796, 0.02279772022797721, 0.03139686031396857, 0.03649635036496359, 0.027497250274972553, 0.029397060293970528, 0.0, 0.04889511048895102, 0.9961003899610038, 0.025497450254974456, 0.023797620237976252, 0.0, 0.03649635036496353, 0.04039596040395971, 0.02769723027697237, 0.03269673032696727, 0.026597340265973404, 0.03239676032396761, 0.03559644035596425, 0.03569643035696427, 0.06939306069393102, 0.032096790320968016, 0.9961003899610041, 0.9963003699630038, 0.03879612038796125, 0.0258974102589741, 0.0, 0.02789721027897221, 0.9966003399660035, 0.03679632036796329, 0.02239776022397762, 0.0, 0.9972002799720028, 0.14988501149885, 0.022097790220977943, 0.9965003499650041, 0.1909809019098084, 0.02669733026697337)
  val dmind = new DenseMatrix(20, 20, ind)

  override def variableSelection(info: InitialInfo): FullStateList = {

    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels) //Not used in SymmetricMain implementation
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels) //Not used in SymmetricMain implementation
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation

    calculateNewState(info.noOfIter, info, FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus), FullStateList(List(FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus))))
  }

  // Update mu and tau
  // helper function for mu tau
  override def nextmutau(oldfullState: FullState, info: InitialInfo): FullState= {
    val prevtau = oldfullState.mt(1)
    val prevmu = oldfullState.mt(0)
    val varMu = 1.0 / (info.tau0 + info.N * prevtau) //the variance for mu
    val meanMu = (info.mu0 * info.tau0 + prevtau * (info.SumObs - sumAllMainInterEff(info.structure, oldfullState.zcoefs, info.zetaLevels, oldfullState.thcoefs, oldfullState.indics))) * varMu
    val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
    val newtau = breeze.stats.distributions.Gamma(info.a + info.N / 2.0, 1.0 / (info.b + 0.5 * YminusMuAndEffects(info.structure, prevmu, oldfullState.zcoefs, oldfullState.thcoefs, oldfullState.indics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
    //oldfullState.copy(mt=DenseVector(newmu,newtau))
    oldfullState.copy(mt=DenseVector(5.3,1.0))
  }

  // Update taus (taua, taub, tauInt)
  // helper function for taus
  override def nexttaus(oldfullState: FullState, info: InitialInfo):FullState= {

    //todo: check if acoef non set values create an issue
    var sumzj = 0.0

    oldfullState.zcoefs.foreachValue( zcoef => {
      sumzj += pow(zcoef - info.alphaPriorMean, 2)
    })

    //todo: check if thcoef non set values create an issue
    var sumThetajk = 0.0
    oldfullState.thcoefs.foreachValue(thcoef => {
      sumThetajk += pow(thcoef -info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })

    val njk = info.structure.sizeOfStructure() // Number of levels of interactions
    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.zetaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    //oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
    oldfullState.copy(tauabth = DenseVector(0.39, 0.22))
  }

  override def nextCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    nextZetaCoefs(oldfullState, info)
  }

  // Update alpha coefficients
  // helper function for alpha coeffs
  def nextZetaCoefs(oldfullState: FullState, info: InitialInfo):FullState={

    val curZetaEstim = DenseVector.zeros[Double](info.zetaLevels)
    (0 until info.zetaLevels).foreach( item => { //For each existing? zeta
      val j = item
      val SXZetaj = info.structure.calcZetaSum(j) // the sum of the observations that have zeta == j on either or both sides
      val Nj = info.structure.calcZetaLength(j) // the number of the observations that have zeta == j on either or both sides
      val SumZeta = sumEffectsOfOtherZetas(info.structure, j, oldfullState.zcoefs) //the sum of the other zeta effects given zeta, for which the given z is on either side (but not on both sides)
      val SinterZeta = sumInterEffGivenZeta(info.structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given zeta, for which the given z is on either side, both sides included
      val varPzeta = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for alphaj
      val meanPzeta = (info.alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXZetaj - Nj * oldfullState.mt(0) - SumZeta - SinterZeta)) * varPzeta //the mean for alphaj
      curZetaEstim(j) = breeze.stats.distributions.Gaussian(meanPzeta, sqrt(varPzeta)).draw()
    })

    oldfullState.copy(zcoefs = curZetaEstim)
    //oldfullState.copy(zcoefs = DenseVector(0.156,0.933,-0.026, 1.357, 2.140, -1.499, 1.270, 1.195, -1.114, -2.619, -1.375, 0.237, 3.077, 0.818, -4.141, 0.567, 1.065, -0.237, -1.294, 1.023 ))
  }

  // Update indicators, interactions and final interaction coefficients
  //Helper function for indicators and interactions

  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo):FullState= {

    val curIndicsEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val curThetaEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    var count = 0.0

    info.structure.foreach( item => {
      val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
      val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b)
      val logInitExp = oldfullState.mt(1) * thcoef * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b) + 0.5 * thcoef))
      val logProb0 = log(1.0 - info.p) //The log of the probability I=0
      val logProb1 = log(info.p) + logInitExp //The log of the probability I=1
      val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
      val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
      val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
      val newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
//      val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

      if (newProb0 < u) {
        //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
        curIndicsEstim(item.a, item.b) = 1.0
        count += 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Njk) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b)))) * varPInter
        curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        curIndicsEstim(item.a,item.b) = 0.0
        curThetaEstim(item.a,item.b) = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
      }
    })

    //oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)

    oldfullState.copy(thcoefs = dm, indics = dmind, finalCoefs = curThetaEstim*:*curIndicsEstim )
  }

  /**
    * Add all the zeta effects for all the other zetas for that specific zeta.
    * e.g. updating z1: (1,1),(1,2),(2,1),(1,3),(1,4),(4,1) => Sum the effects for: z2*NoOfObs for that category + z2*NoOfObs for that category + z3*NoOfObs for that category + z4*NoOfObs for that category + z4*NoOfObs for that category
    */
  def sumEffectsOfOtherZetas(structure: DVStructure, zetaIndex: Int, zetaEff: DenseVector[Double]): Double = {
    //returns the element which is not zetaIndex. It doesn't take into account the cases where both sides are zetaIndex because getAllOtherZetasItemsForGivenZ works on a structure that does not involve the (j,j) cases
    def notZeta(k1: Int, k2: Int): Int={
      if(k1!=zetaIndex) k1
      else k2
    }
    structure.getAllOtherZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length * zetaEff(notZeta(elem._1._1, elem._1._2))).reduce(_+_)
  }

  /**
    * Calculate the sum of all the zeta 1 and all the zeta 2 effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, zetaEff: DenseVector[Double], nz: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var suma = 0.0
    var sumb = 0.0
    var sumInter = 0.0

    // For zeta effects in first column
    structure.foreach(item => {
      suma += item.list.length * zetaEff(item.a)
    })

    // For zeta effects in second column
    structure.foreach(item => {
      sumb += item.list.length * zetaEff(item.b)
    })

    // Add all the interaction effects for a given alpha and a given beta taking advantage of the DVStructure
    structure.foreach( item => {
      sumInter += item.list.length * indics(item.a, item.b) * interEff(item.a, item.b)
    })

    suma + sumb + sumInter
  }

  /**
    * Add all the interaction effects for a given zeta. Adds all the interactions for which zeta is on either side. Includes the doubles bcs getZetasItemsForGivenZ uses a structure that includes everything
    */
  def sumInterEffGivenZeta(structure: DVStructure, zetaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    structure.getZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length * indics(elem._1._1, elem._1._2) * interEff(elem._1._1, elem._1._2)).reduce(_+_)
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(structure:DVStructure, mu: Double, zetaEff: DenseVector[Double], interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sum = 0.0

    structure.foreach( item => {
      val a = item.a
      val b = item.b
      sum += item.list.map(x => scala.math.pow(x - mu - zetaEff(a) - zetaEff(b) - interEff(a, b) * indics(a, b), 2)).sum
    })
    sum
  }

  override def printResults(statesResults: FullStateList): Unit = {
    println("zetas")
    val zcoefficients = statesResults.fstateL.map(f => f.zcoefs)
    //println(acoefficients)
    val zcoefMat = DenseMatrix(zcoefficients.map(_.toArray): _*)
    val meanValsZcoef = mean(zcoefMat(::, *))
    println(meanValsZcoef)

    val matrices = calculateAndPrintCommons(statesResults)

    // Save the results to a csv file
    val noOfInters = matrices(2).cols
    val mergedMatrix = DenseMatrix.horzcat(matrices(0), matrices(1), zcoefMat, matrices(2), matrices(3))
    saveToCSV(mergedMatrix, getFileNameToSaveResults("allcoefs"))
    //    saveToCSV(matrices(0), getFileNameToSaveResults("mutau"))
    //    saveToCSV(matrices(1), getFileNameToSaveResults("taus"))
    //    saveToCSV(zcoefMat, getFileNameToSaveResults("zetas"))
    //    saveToCSV(matrices(2), getFileNameToSaveResults("thetas"))
    //    saveToCSV(matrices(2)(::, 0 to noOfInters/2), getFileNameToSaveResults("thetas1"))
    //    saveToCSV(matrices(2)(::, (noOfInters/2)+1 to noOfInters-1), getFileNameToSaveResults("thetas2"))
    //    saveToCSV(matrices(3)(::, 0 to noOfInters/2), getFileNameToSaveResults("indics1"))
    //    saveToCSV(matrices(3)(::, (noOfInters/2)+1 to noOfInters-1), getFileNameToSaveResults("indics2"))
    //    saveToCSV(matrices(3), getFileNameToSaveResults("indics"))
  }

  override protected def getFileNameToSaveResults(param: String): String = {
    val filePath = getMainFilePath.concat("/symmetricMain3-10mScalaRestryAllConstFindZ-")
    val pathToFiles = Map("mutau" -> filePath.concat("mutau.csv") ,
      "taus" -> filePath.concat("taus.csv"),
      "zetas" -> filePath.concat("zetas.csv"),
      "thetas" -> filePath.concat("thetas.csv"),
      "indics" -> filePath.concat("indics.csv"),
      "allcoefs" -> filePath.concat("allcoefs.csv")
    )
    pathToFiles(param)
  }

  override def getInputFilePath(): String = getMainFilePath.concat("/simulInterSymmetricMain.csv")

}