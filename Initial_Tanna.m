% Load the measured source and its sample rate
load ('~/Desktop/reverberation/room-measures/source.mat')
% Load the measured IR
addr1a = '~/Desktop/reverberation/room-measures/tanna/tanna1a.mat';
addr1b = '~/Desktop/reverberation/room-measures/tanna/tanna1b.mat';
addr2a = '~/Desktop/reverberation/room-measures/tanna/tanna2a.mat';
addr2b = '~/Desktop/reverberation/room-measures/tanna/tanna2b.mat';
addr3a = '~/Desktop/reverberation/room-measures/tanna/tanna3a.mat';
addr3b = '~/Desktop/reverberation/room-measures/tanna/tanna3b.mat';
addr4 = '~/Desktop/reverberation/room-measures/tanna/tanna4.mat';
addr5 = '~/Desktop/reverberation/room-measures/tanna/tanna5.mat';
addr6a = '~/Desktop/reverberation/room-measures/tanna/tanna6a.mat';
addr6b = '~/Desktop/reverberation/room-measures/tanna/tanna6b.mat';
addr6c = '~/Desktop/reverberation/room-measures/tanna/tanna6c.mat';

% Space volume
% V = 1.6187e+03-70.4700-166.0500-70.2;
V = 1.3120e+03;

% Different source and receiver positions in the model
dLaser = 0.11;
H = 7.11;
src =[6  15.91-dLaser  2.37; ...
    6  15.91-dLaser  2.37; ...
    6  15.91-dLaser  2.37; ...
    ];
rcv = [
    6.4+dLaser    10.9-dLaser  2.21; ...
    4.295+dLaser  7.375-dLaser 3.41; ...
    7.54+dLaser   10.03-dLaser 2.51; ...
    ];
ref_direc = [
    0   -1   0; ...
    ];

% Wall absorption coefficients at different frequency bands
f = [125 250 500 1000 2000 4000 8000];
Gypsum_board = [0.3 0.69 1 0.81 0.66 0.62 0.62];
wood = [0.09 0.06 0.05 0.05 0.05 0.04 0.04];
StageWood = [0.1 0.07 0.06 0.06 0.06 0.06 0.06];
steel_frame = [0.15 0.1 0.06 0.04 0.04 0.05 0.05];
Acoustical_plaster = [0.17 0.36 0.66 0.65 0.62 0.68 0.68];
RPG_QRD = [0.06 0.15 0.45 0.95 0.88 0.91 0.91];
chairs = [0.49 0.66 0.8 0.88 0.82 0.70 0.70];

% f_Gypsum_board = Wall_Filter(f, Gypsum_board, fs,  DispFilter);
% f_wood = Wall_Filter(f, wood, fs,  DispFilter);
% f_StageWood = Wall_Filter(f, StageWood, fs, DispFilter);
% f_steel_frame = Wall_Filter(f, steel_frame, fs, DispFilter);
% f_Acoustical_plaster = Wall_Filter(f, Acoustical_plaster, fs,  DispFilter);
% f_RPG_QRD = Wall_Filter(f, RPG_QRD, fs,  DispFilter);

f_Gypsum_board = sqrt(1-Gypsum_board(Band_num));
f_StageWood = sqrt(1-StageWood(Band_num));
f_steel_frame = sqrt(1-steel_frame(Band_num));
f_Acoustical_plaster = sqrt(1-Acoustical_plaster(Band_num));
f_RPG_QRD = sqrt(1-RPG_QRD(Band_num));
f_chairs = sqrt(1-chairs(Band_num));

beta = [f_Gypsum_board;...     % 1
    f_Gypsum_board;...
    f_Gypsum_board;...     
    f_StageWood;...
    f_StageWood;...       5
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...
    f_Gypsum_board;...
    f_Gypsum_board;...       10
    f_Gypsum_board;...
    f_RPG_QRD;...
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...       15
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...       20
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...
    f_StageWood;...       25
    f_StageWood;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       30
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       35
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       40
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       45
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       50
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       55
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       60
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       65
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       70
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_chairs;...       75
    f_chairs;...
    f_chairs;...
    f_chairs;...
    f_steel_frame;...
    f_steel_frame;...       80
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...       85
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...       90
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_steel_frame;...
    f_Acoustical_plaster;...       95
    f_RPG_QRD;...
    f_RPG_QRD;...
    f_RPG_QRD;...
    f_RPG_QRD;...
    f_RPG_QRD;...       100
    f_StageWood;...
    f_Gypsum_board;...
    f_RPG_QRD;...
    f_RPG_QRD;...
    f_RPG_QRD;...       105
    f_RPG_QRD;...
    f_RPG_QRD;...
    ];

vertex = [
    0     14.4    H; ...
    1     20      H; ...
    11    20      H; ...
    12    14.4    H; ...
    0     14.4    0; ...        5
    0     15.2    0; ...
    0     15.2    0.72; ...
    1     20      0.72; ...
    11    20      0.72; ...
    12    15.2    0.72; ...     10
    12    15.2    0; ...
    12    14.4    0; ...
    1.5   15.2    0.72; ...
    1.5   14.4    0.72; ...
    10.5  14.4    0.72; ...     15
    10.5  15.2    0.72; ...
    1.5   15.2    0; ...
    1.5   14.4    0; ...
    10.5  14.4    0; ...
    10.5  15.2    0; ...        20
    0     0       H; ...
    12    0       H; ...
    0     0       0; ...
    12    0       0; ...
    % stage stairs
    0     14.4    0; ...        25
    1.5   14.4    0; ...
    0     14.4    0.18; ...
    1.5   14.4    0.18; ...
    0     14.65   0.18; ...
    1.5   14.65   0.18; ...      30
    0     14.65   0.36; ...
    1.5   14.65   0.36; ...
    0     14.9    0.36; ...
    1.5   14.9    0.36; ...
    0     14.9    0.54; ...      35
    1.5   14.9    0.54; ...
    0     15.2   0.54; ...
    1.5   15.2   0.54; ...
    10.5  14.4    0; ...
    12    14.4    0; ...        40
    10.5  14.4    0.18; ...
    12    14.4    0.18; ...
    10.5  14.65   0.18; ...
    12    14.65   0.18; ...
    10.5  14.65   0.36; ...     45
    12    14.65   0.36; ...
    10.5  14.9    0.36; ...
    12    14.9    0.36; ...
    10.5  14.9    0.54; ...
    12    14.9    0.54; ...     50
    10.5  15.2   0.54; ...
    12    15.2   0.54; ...
    % stairs
    1.65  12.75   0; ...
    10.35 12.75   0; ...
    1.65  12.75   0.15; ...     55
    10.35 12.75   0.15; ...
    1.65  12.3    0.15; ...
    10.35 12.3    0.15; ...
    1.65  12.3    0.3; ...
    10.35 12.3    0.3; ...      60
    1.65  11.85   0.3; ...
    10.35 11.85   0.3; ...
    1.65  11.85   0.45; ...
    10.35 11.85   0.45; ...
    1.65  11.4    0.45; ...     65
    10.35 11.4    0.45; ...
    1.65  11.4    0.6; ...
    10.35 11.4    0.6; ...
    1.65  10.95   0.6; ...
    10.35 10.95   0.6; ...      70
    1.65  10.95   0.75; ...
    10.35 10.95   0.75; ...
    1.65  10.5    0.75; ...
    10.35 10.5    0.75; ...
    1.65  10.5    0.9; ...      75
    10.35 10.5    0.9; ...
    1.65  10.05   0.9; ...
    10.35 10.05   0.9; ...
    1.65  10.05   1.05; ...
    10.35 10.05   1.05; ...     80
    1.65  9.6     1.05; ...
    10.35 9.6     1.05; ...
    1.65  9.6     1.2; ...
    10.35 9.6     1.2; ...
    1.65  9.15    1.2; ...      85
    10.35 9.15    1.2; ...
    1.65  9.15    1.35; ...
    10.35 9.15    1.35; ...
    1.65  8.7     1.35; ...
    10.35 8.7     1.35; ...     90
    1.65  8.7     1.5; ...
    10.35 8.7     1.5; ...
    1.65  8.25    1.5; ...
    10.35 8.25    1.5; ...
    1.65  8.25    1.65; ...     95
    10.35 8.25    1.65; ...
    1.65  7.8     1.65; ...
    10.35 7.8     1.65; ...
    1.65  7.8     1.8; ...
    10.35 7.8     1.8; ...      100
    1.65  7.35    1.8; ...
    10.35 7.35    1.8; ...
    1.65  7.35    1.95; ...
    10.35 7.35    1.95; ...
    1.65  6.9     1.95; ...     105
    10.35 6.9     1.95; ...
    1.65  6.9     2.1; ...
    10.35 6.9     2.1; ...
    1.65  6.45    2.1; ...
    10.35 6.45    2.1; ...      110
    1.65  6.45    2.25; ...
    10.35 6.45    2.25; ...
    1.65  6       2.25; ...
    10.35 6       2.25; ...
    1.65  6       2.4; ...      115
    10.35 6       2.4; ...
    
    0     6       2.4; ...
    12    6       2.4; ...
    0     5.55    2.4; ...
    12    5.55    2.4; ...      120
    
    0     5.55    2.55; ...
    12    5.55    2.55; ...
    0     5.1     2.55; ...
    12    5.1     2.55; ...
    0     5.1     2.7; ...      125
    12    5.1     2.7; ...
    0     4.65    2.7; ...
    12    4.65    2.7; ...
    0     4.65    2.85; ...
    12    4.65    2.85; ...     130
    0     4.2     2.85; ...
    12    4.2     2.85; ...
    0     4.2     3; ...
    12    4.2     3; ...
    0     3.75    3; ...        135
    12    3.75    3; ...
    0     3.75    3.15; ...
    12    3.75    3.15; ...
    0     3.3     3.15; ...
    12    3.3     3.15; ...     140
    0     3.3     3.3; ...
    12    3.3     3.3; ...
    0     2.85    3.3; ...
    12    2.85    3.3; ...
    0     2.85    3.45; ...     145
    12    2.85    3.45; ...
    0     2.4     3.45; ...
    12    2.4     3.45; ...
    0     2.4     3.6; ...
    12    2.4     3.6; ...      150
    0     1.95    3.6; ...
    12    1.95    3.6; ...
    0     1.95    3.75; ...
    12    1.95    3.75; ...
    0     1.5     3.75; ...     155
    12    1.5     3.75; ...
    0     1.5     3.9; ...
    12    1.5     3.9; ...
    0     0       3.9; ...
    12    0       3.9; ...      160
    % stair handrail left
    1.65  13.2    0; ...
    1.65  13.2    1.07; ...
    1.65  12.88   1.07; ...
    1.65  6.5     3.17; ...
    1.65  6       3.17; ...     165
    1.65  6       0; ...
    0     6       0; ...
    0     6       3.17; ...
    1.6   13.2    0; ...
    1.6   13.2    1.07; ...     170
    1.6   12.88   1.07; ...
    1.6   6.5     3.17; ...
    1.6   6.05    3.17; ...
    1.6   6.05    0; ...
    0     6.05    0; ...        175
    0     6.05    3.17; ...
    % stair handrail right
    10.35 13.2    0; ...
    10.35 13.2    1.07; ...
    10.35 12.88   1.07; ...
    10.35 6.5     3.17; ...     180
    10.35 6       3.17; ...
    10.35 6       0; ...
    12    6       0; ...
    12    6       3.17; ...
    10.4  13.2    0; ...        185
    10.4  13.2    1.07; ...
    10.4  12.88   1.07; ...
    10.4  6.5     3.17; ...
    10.4  6.05    3.17; ...
    10.4  6.05    0; ...        190
    12    6.05    0; ...
    12    6.05    3.17; ...
    % wooden wall
    0.05    0      2.05; ...
    0.05    14.4   2.05; ...
    1.05    19.95  2.05; ...        195
    10.95   19.95  2.05; ...
    11.95   14.4   2.05; ...
    11.95   0      2.05; ...
    0.05    0      H; ...
    0.05    14.4   H; ...           200
    1.05    19.95  H; ...
    10.95   19.95  H; ...
    11.95   14.4   H; ...
    11.95   0      H; ...
    % added points
    1     20      0; ...            205
    11    20      0; ...
    0     0      2.05; ...
    0     14.4   2.05; ...
    1     20     2.05; ...     
    11    20     2.05; ...        210
    12    14.4   2.05; ...
    12    0      2.05; ...
    ];

wall = [
    1  2  205 5  0   0   0   0; ... % 1
    2  3  9  8   0   0   0   0; ... % 2
    3  4  12 206 0   0   0   0; ... % 3
    7  8  9  10  0   0   0   0; ... % 4
    7  13 17 6   0   0   0   0; ... % 5
    16 10 11 20  0   0   0   0; ... % 6
    13 14 18 17  0   0   0   0; ... % 7
    15 16 20 19  0   0   0   0; ... % 8
    5  12 24 23  0   0   0   0; ... % 9
    1  5  23 21  0   0   0   0; ... % 10
    4  22 24 12  0   0   0   0; ...
    21 23 24 22  0   0   0   0; ...
    % stage stairs
    25 27 28 26  0   0   0   0; ...
    27 29 30 28  0   0   0   0; ...
    29 31 32 30  0   0   0   0; ... % 15
    31 33 34 32  0   0   0   0; ...
    33 35 36 34  0   0   0   0; ...
    35 37 38 36  0   0   0   0; ...
    37 7  13 38  0   0   0   0; ...
    39 41 42 40  0   0   0   0; ... % 20
    41 43 44 42  0   0   0   0; ...
    43 45 46 44  0   0   0   0; ...
    45 47 48 46  0   0   0   0; ...
    47 49 50 48  0   0   0   0; ...
    49 51 52 50  0   0   0   0; ... % 25
    51 16 10 52  0   0   0   0; ...
    % stairs audience
    53 54 56 55  0   0   0   0; ...
    55 56 58 57  0   0   0   0; ...
    57 58 60 59  0   0   0   0; ...
    59 60 62 61  0   0   0   0; ... % 30
    61 62 64 63  0   0   0   0; ...
    63 64 66 65  0   0   0   0; ...
    65 66 68 67  0   0   0   0; ...
    67 68 70 69  0   0   0   0; ...
    69 70 72 71  0   0   0   0; ... % 35
    71 72 74 73  0   0   0   0; ...
    73 74 76 75  0   0   0   0; ...
    75 76 78 77  0   0   0   0; ...
    77 78 80 79  0   0   0   0; ...
    79 80 82 81  0   0   0   0; ... % 40
    81 82 84 83  0   0   0   0; ...
    83 84 86 85  0   0   0   0; ...
    85 86 88 87  0   0   0   0; ...
    87 88 90 89  0   0   0   0; ...
    89 90 92 91  0   0   0   0; ... % 45
    91 92 94 93  0   0   0   0; ...
    93 94 96 95  0   0   0   0; ...
    95 96 98 97  0   0   0   0; ...
    97 98 100 99 0   0   0   0; ...
    99 100 102 101 0 0   0   0; ... % 50
    101 102 104 103 0 0  0   0; ...
    103 104 106 105 0 0  0   0; ...
    105 106 108 107 0 0  0   0; ...
    107 108 110 109 0 0  0   0; ...
    109 110 112 111 0 0  0   0; ... % 55
    111 112 114 113 0 0  0   0; ...
    113 114 116 115 0 0  0   0; ...
    117 118 120 119 0 0  0   0; ...
    119 120 122 121 0 0  0   0; ...
    121 122 124 123 0 0  0   0; ... % 60
    123 124 126 125 0 0  0   0; ...
    125 126 128 127 0 0  0   0; ...
    127 128 130 129 0 0  0   0; ...
    129 130 132 131 0 0  0   0; ...
    131 132 134 133 0 0  0   0; ... % 65
    133 134 136 135 0 0  0   0; ...
    135 136 138 137 0 0  0   0; ...
    137 138 140 139 0 0  0   0; ...
    139 140 142 141 0 0  0   0; ...
    141 142 144 143 0 0  0   0; ... % 70
    143 144 146 145 0 0  0   0; ...
    145 146 148 147 0 0  0   0; ...
    147 148 150 149 0 0  0   0; ...
    149 150 152 151 0 0  0   0; ...
    151 152 154 153 0 0  0   0; ... % 75
    153 154 156 155 0 0  0   0; ...
    155 156 158 157 0 0  0   0; ...
    157 158 160 159 0 0  0   0; ...
    % stair handrail left
    161 166 165 164 163 162 0 0; ...
    165 166 167 168 0   0   0 0; ...% 80
    169 170 171 172 173 174 0 0; ...
    176 175 174 173 0   0   0 0; ...
    161 162 170 169 0   0   0 0; ...
    162 163 171 170 0   0   0 0; ...
    163 164 172 171 0   0   0 0; ...% 85
    164 165 168 176 173 172 0 0; ...
    % stair handrail right
    177 178 179 180 181 182 0 0; ...
    184 183 182 181 0   0   0 0; ...
    185 190 189 188 187 186 0 0; ...
    189 190 191 192 0   0   0 0; ...% 90
    185 186 178 177 0   0   0 0; ...
    186 187 179 178 0   0   0 0; ...
    187 188 180 179 0   0   0 0; ...
    180 188 189 192 184 181 0 0; ... 94
    % ceiling
    1   21  22  4   3   2   0 0; ... 95
    % wooden walls
    199 200 194 193 0   0   0 0; ...
    200 201 195 194 0   0   0 0; ...
    201 202 196 195 0   0   0 0; ...
    202 203 197 196 0   0   0 0; ...
    203 204 198 197 0   0   0 0; ... 100
    % additional stage
    16  15  14  13  0   0   0  0; ...
    14  15  19  18  0   0   0  0; ...
    193 194 208 207 0   0   0  0; ...
    194 195 209 208 0   0   0  0; ...
    195 196 210 209 0   0   0  0; ... 105
    196 197 211 210 0   0   0  0; ...
    197 198 212 211 0   0   0  0; ...
    ];

wnum = size(wall,1); % calculate the total wall number


