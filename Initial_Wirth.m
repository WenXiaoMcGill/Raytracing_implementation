% Load the measured source and its sample rate
load ('~/Desktop/reverberation/room-measures/source.mat')
% Load the measured IR
addr1a = '~/Desktop/reverberation/room-measures/wirth/Wirth1a.mat';
addr1b = '~/Desktop/reverberation/room-measures/wirth/Wirth1b.mat';
addr2a = '~/Desktop/reverberation/room-measures/wirth/Wirth2a.mat';
addr2b = '~/Desktop/reverberation/room-measures/wirth/Wirth2b.mat';
addr3a = '~/Desktop/reverberation/room-measures/wirth/Wirth3a.mat';
addr3b = '~/Desktop/reverberation/room-measures/wirth/Wirth3b.mat';

% Space volume
% V = 1.5013e+03-4.9603-100.2506;
V = 1.3961e+03;

% Different source and receiver positions in the model
dLaser = 0.11;
H = 5.99;
src =[5.87+dLaser  4.93+dLaser+1.82  1.64; ...
      5.87+dLaser  4.93+dLaser+1.82  1.64; ...
      5.87+dLaser  4.93+dLaser+1.82  1.64; ...
     ];
rcv = [
    5.98+dLaser  4.93+4.76+1.82+dLaser*2  1.65; ...
    9.78+dLaser  12.415+dLaser            1.65; ...
    2.62+dLaser  10.1+dLaser              1.65; ...
    ];
ref_direc = [
    0   1   0; ...
    1   0   0; ...
    -1  0   0; ...
    ];

% Wall absorption coefficients at different frequency bands
f = [125 250 500 1000 2000 4000 8000];
Gypsum_board = [0.3 0.69 1 0.81 0.66 0.62 0.62];
curtain = [0.36 0.26 0.51 0.45 0.62 0.76 0.76];
plaster = [0.02 0.02 0.03 0.03 0.04 0.05 0.05];
StageWood = [0.1 0.07 0.06 0.06 0.06 0.06 0.06];
chair = [0.40 0.50 0.58 0.61 0.58 0.5 0.5];

% f_Gypsum_board = Wall_Filter(f, Gypsum_board, fs,  DispFilter);
% f_curtain = Wall_Filter(f, curtain, fs,  DispFilter);
% f_plaster = Wall_Filter(f, plaster, fs, DispFilter);
% f_StageWood = Wall_Filter(f, StageWood, fs, DispFilter);
% f_chair = Wall_Filter(f, chair, fs, DispFilter);

f_Gypsum_board = sqrt(1-Gypsum_board(Band_num));
f_curtain = sqrt(1-curtain(Band_num));
f_plaster = sqrt(1-plaster(Band_num));
f_StageWood = sqrt(1-StageWood(Band_num));
f_chair = sqrt(1-chair(Band_num));

beta = [
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_curtain;...
    f_plaster;...
    f_StageWood;...       5
    f_plaster;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...       10
    f_plaster;...
    f_plaster;...
    f_plaster;...
    f_plaster;...
    f_plaster;...       15
    f_plaster;...
    f_curtain;...
    f_curtain;...
    f_curtain;...
    f_curtain;...       20
    f_curtain;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...       25
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...
    f_Gypsum_board;...      30
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...             35
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...             40
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...             45
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...             50
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...
    f_chair;...             55
];

vertex = [
    0       0      0; ...
    12.42   0      0; ...
    12.42   20.18  0; ...
    0       20.18  0; ...
    0       0      H; ...       5
    12.42   0      H; ...
    12.42   20.18  H; ...
    0       20.18  H; ...
    % column
    5.61    0.91   0; ...
    6.52    0.91   0; ...       10
    6.52    1.82   0; ...
    5.61    1.82   0; ...
    5.61    0.91   H; ...
    6.52    0.91   H; ...
    6.52    1.82   H; ...       15
    5.61    1.82   H; ...
    % 
    0.93    0.52   3.66; ...    
    11.49   0.52   3.66; ...
    11.49   20.18  3.66; ...
    0.93    20.18  3.66; ...    20
    0.93    0.52   H; ...
    11.49   0.52   H; ...
    11.49   20.18  H; ...
    0.93    20.18  H; ...
    0       0      3.66; ...    25
    0.93    0      3.66; ...
    0.93    20.18  3.66; ...
    0       20.18  3.66; ...
    11.49   0      3.66; ...
    12.42   0      3.66; ...    30
    12.42   20.18  3.66; ...
    11.49   20.18  3.66; ...
    % curtain
    1.21    2.5    2.62; ...
    1.21    16.5   2.62; ...
    1.21    16.5   H; ...       35
    1.21    2.5    H; ...
    0.93    2.5    2.62; ...
    0.93    16.5   2.62; ...
    0.93    16.5   H; ...     
    0.93    2.5    H; ...       40
    % stairs
    2.87    16.53  0.2; ...
    9.7     16.53  0.2; ...
    2.87    17.75  0.2; ...
    9.7     17.75  0.2; ...  
    2.87    16.53  0; ...       45
    9.7     16.53  0; ...
    2.87    17.75  0; ...
    9.7     17.75  0; ...
    
    2.87    17.75  0.4; ...
    9.7     17.75  0.4; ...     50
    2.87    18.96  0.4; ...
    9.7     18.96  0.4; ...
    2.87    17.75  0; ...
    9.7     17.75  0; ...
    2.87    18.96  0; ...       55
    9.7     18.96  0; ...
    
    2.87    18.96  0.6; ... 
    9.7     18.96  0.6; ...
    2.87    20.18  0.6; ...
    9.7     20.18  0.6; ...     60
    2.87    18.96  0; ... 
    9.7     18.96  0; ...
    2.87    20.18  0; ...
    9.7     20.18  0; ...
    
    % chairs
    2.87    15.31  0; ...       65
    9.7     15.31  0; ...     
    2.87    15.31  0.92; ...
    9.7     15.31  0.92; ...
    2.87    15.01  0; ... 
    9.7     15.01  0; ...       70
    2.87    15.01  0.92; ...
    9.7     15.01  0.92; ...
    
    2.87    16.53  0; ...
    9.7     16.53  0; ...  
    2.87    16.53  0.92; ...    75
    9.7     16.53  0.92; ...
    2.87    16.23  0; ...
    9.7     16.23  0; ...  
    2.87    16.23  0.92; ...  
    9.7     16.23  0.92; ...    80
    
    2.87    17.75  0.2; ...
    9.7     17.75  0.2; ...
    2.87    17.75  1.12; ...
    9.7     17.75  1.12; ...  
    2.87    17.45  0.2; ...     85
    9.7     17.45  0.2; ...
    2.87    17.45  1.12; ...
    9.7     17.45  1.12; ...  
    
    2.87    18.96  0.4; ...
    9.7     18.96  0.4; ...     90
    2.87    18.96  1.32; ...
    9.7     18.96  1.32; ... 
    2.87    18.66  0.4; ...
    9.7     18.66  0.4; ...
    2.87    18.66  1.32; ...    95
    9.7     18.66  1.32; ... 
    
    2.87    20.18  0.6; ...  
    9.7     20.18  0.6; ...
    2.87    20.18  1.52; ...
    9.7     20.18  1.52; ...    100
    2.87    19.88  0.6; ...  
    9.7     19.88  0.6; ...
    2.87    19.88  1.52; ...
    9.7     19.88  1.52; ...
    ];

wall = [
    1   2   6   5; ...
    2   3   7   6; ...
    3   4   8   7; ...
    4   1   5   8; ...
    4   3   2   1; ...          5
    5   6   7   8; ...
    % column
    13  14  10  9;...
    14  15  11  10;...
    15  16  12  11;...
    16  13  9   12;...         10
    %
    17  18  22  21; ...
    18  19  23  22; ...
    20  17  21  24; ...
    25  26  27  28; ...
    26  29  18  17; ...         15
    29  30  31  32; ...
    % curtain
    36  35  34  33; ...
    37  38  39  40; ...
    33  37  40  36; ...
    33  34  38  37; ...         20
    34  35  39  38; ...
    % stairs
    43  44  42  41; ...
    42  44  48  46; ...
    43  41  45  47; ...         
    51  52  50  49; ...         25
    50  52  56  54; ...
    51  49  53  55; ...
    59  60  58  57; ...
    58  60  64  62; ...         
    59  57  61  63; ...         30
    % chairs
    68  67  65  66; ...
    71  72  70  69; ...
    67  71  69  65; ...
    67  68  72  71; ...
    70  72  68  66; ...         35
    
    76  75  73  74; ...
    79  80  78  77; ...
    75  79  77  73; ...
    75  76  80  79; ...
    78  80  76  74; ...         40
    
    84  83  81  82; ...
    87  88  86  85; ...
    83  87  85  81; ...
    83  84  88  87; ...
    86  88  84  82; ...         45
    
    92  91  89  90; ...
    95  96  94  93; ...
    91  95  93  89; ...
    91  92  96  95; ...
    94  96  92  90; ...         50
    
    100  99  97  98; ...
    103  104 102 101; ...
    99  103  101  97; ...
    99  100  104  103; ...
    102 104  100  98; ...         55
];

wnum = size(wall,1); % calculate the total wall number
