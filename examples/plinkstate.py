import math

class pniniVars :
    
    #### DEFAULTS:
    DEFAULT_KUHN    = 0.88
    DEFAULT_CONTOUR = 0.38
    DEFAULT_MAX_AA  = 200
    
    CELLULOSE_STEP_X   = 1.026
    MYOSIN_STEP_X      = 8.0
    MYOSIN_HALF_STEP_X = 4.0
    CELL_STEP2_X       = CELLULOSE_STEP_X * 2.0 # 2.052
    CELL_STEP3_X       = CELLULOSE_STEP_X * 3.0 # 3.078
    CELL_STEP4_X       = CELLULOSE_STEP_X * 4.0 # 4.104
    CELL_STEP5_X       = CELLULOSE_STEP_X * 5.0 # 5.130
    CELL_STEP6_X       = CELLULOSE_STEP_X * 6.0 # 6.156
    CELL_STEP7_X       = CELLULOSE_STEP_X * 7.0 # 7.182
    CELL_STEP8_X       = CELLULOSE_STEP_X * 8.0 # 8.208
    
    CELLULOSE_STEP_Y = 0.0
    CELLULOSE_SHEET_STEP_Y = 2.0
    MYOSIN_STEP_Y    = 0.0
    
    IS_SHEET_SUBSTRATE = False
    
    CEL7A_KON      = 31.0
    CEL7A_KON_M    = 1.0e-3
    KIN_KON        = 20.0
    KIN_KON_M      = 1e-6
    
    DEFAULT_MIN_DOCKED     = 1
    DEFAULT_MAX_DOCKED     = 8
    DEFAULT_DOCK_STEP      = 1
    DEFAULT_DOCKPERCENT    = 1.0
    DOCKPERCENT_LOW        = 0.4
    DOCKPERCENT_MID        = 0.8
    DOCKPERCENT_HIGH       = 1.0    
    
    DOCKLIST_DEFAULT = []
    DOCKLIST_NODOCK  = [0]
    DOCKLIST_HIGH    = [0,3,5,7]
    DOCKLIST_HIGH2   = [0,2,4,6]
    
    DEFAULT_MIN_SEGMENT_NUM = 10
    DEFAULT_MAX_SEGMENT_NUM = 121
    DEFAULT_SEGMENT_STEP    = 1
        
    PRESETGEOM_CELLULASE = 'cellulase'
    PRESETGEOM_NULL      = 'null'
    DEFAULT_PRESET_GEOM  = 'cellulase'
    
    YLIM_MIN_DEFAULT = None
    YLIM_MAX_DEFAULT = None
    
    XLIM_MIN_DEFAULT = None
    XLIM_MAX_DEFAULT = None
    
    ############################
    Kcat_TIMES = []
    Lower_LL   = []
    Higher_LL  = []
    Kcat_TIMES_LOWER = []
    
    ############################
    
    KUHN    = None
    CONTOUR = DEFAULT_CONTOUR
    
    USE_AA  = True
    SET_AA  = None
    MAX_AA  = DEFAULT_MAX_AA
    
    etol_abs = None
    etol_rel = None
    limit_num = None
    cycle_num = None
    option_dict = None
    stepXNum = None
    stepXHNum = None
    stepYZNum = None
    stepYZHNum = None
    
    STEP_X = None
    STEP_Y = None
    
    MIN_TARGET_STEP_X = None
    MAX_TARGET_STEP_X = None
    MIN_BACK_STEP_X   = None
    MAX_BACK_STEP_X   = None
    MIN_TARGET_STEP_Y = None
    MAX_TARGET_STEP_Y = None
    MIN_BACK_STEP_Y   = None
    MAX_BACK_STEP_Y   = None
    
    SHEET_TARGET_DEGREE  = None
    SHEET_NEUTRAL_DEGREE = None
    SHEET_BACK_DEGREE    = None
    
    KON = None
    KON_M = None
    
    MIN_DOCKED = None
    MAX_DOCKED = None
    DOCK_STEP = None
    DOCK_STEP_NUM = None
    DOCKPERCENT = DEFAULT_DOCKPERCENT

    DOCKLIST = None
    
    MIN_SEGMENT_NUM = None
    MAX_SEGMENT_NUM = None
    SEGMENT_STEP = None
    
    PRESETGEOM = None
    
    MARKSIZE = None

    YLIM_MIN = None
    YLIM_MAX = None
    
    XLIM_MIN = None
    XLIM_MAX = None
    
    HLINE           = 3e-3
    HLINE_THICKNESS = 2.0
    DRAW_HLINE      = True
    
    IMAGE_PATH      = 'img/'
    IMAGE_NAME      = ''
    IMAGE_JPG       = ''
    IMAGE_FULL_JPG  = None
    IMAGE_FULL_NAME = None
    
    PRESET_001 = '001_geom_null_cellstep_linear_kcel_nodock_kuhn088_target1_noback'
    PRESET_002 = '002_geom_null_cellstep_linear_kcel_nodock_kuhn088_target2_noback'
    PRESET_003 = '003_geom_null_cellstep_linear_kcel_nodock_kuhn088_target3_noback'
    PRESET_004 = '004_geom_null_cellstep_linear_kcel_nodock_kuhn088_target4_noback'
    PRESET_005 = '005_geom_null_cellstep_linear_kcel_nodock_kuhn088_target5_noback'
    PRESET_006 = '006_geom_null_cellstep_linear_kcel_nodock_kuhn088_target6_noback'
    
    PRESET_007 = '007_geom_null_cellstep_linear_kcel_nodock_kuhn176_target1_noback'
    PRESET_008 = '008_geom_null_cellstep_linear_kcel_nodock_kuhn176_target2_noback'
    PRESET_009 = '009_geom_null_cellstep_linear_kcel_nodock_kuhn176_target3_noback'
    PRESET_010 = '010_geom_null_cellstep_linear_kcel_nodock_kuhn176_target4_noback'
    PRESET_011 = '011_geom_null_cellstep_linear_kcel_nodock_kuhn176_target5_noback'
    PRESET_012 = '012_geom_null_cellstep_linear_kcel_nodock_kuhn176_target6_noback'
    
    PRESET_013 = '013_geom_null_cellstep_linear_kcel_nodock_kuhn088_target1on_noback'
    PRESET_014 = '014_geom_null_cellstep_linear_kcel_nodock_kuhn088_target2on_noback'
    PRESET_015 = '015_geom_null_cellstep_linear_kcel_nodock_kuhn088_target3on_noback'
    PRESET_016 = '016_geom_null_cellstep_linear_kcel_nodock_kuhn088_target4on_noback'
    PRESET_017 = '017_geom_null_cellstep_linear_kcel_nodock_kuhn088_target5on_noback'
    PRESET_018 = '018_geom_null_cellstep_linear_kcel_nodock_kuhn088_target6on_noback'
    
    PRESET_019 = '019_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target1_noback'
    PRESET_020 = '020_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target2_noback'
    PRESET_021 = '021_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target3_noback'
    PRESET_022 = '022_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target4_noback'
    PRESET_023 = '023_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target5_noback'
    PRESET_024 = '024_geom_cell_cellstep_linear_kcel_nodock_kuhn088_target6_noback'
    
    PRESET_025 = '025_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target1_noback'
    PRESET_026 = '026_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target2_noback'
    PRESET_027 = '027_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target3_noback'
    PRESET_028 = '028_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target4_noback'
    PRESET_029 = '029_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target5_noback'
    PRESET_030 = '030_geom_cell_cellstep_linear_kcel_dock0456_kuhn088_target6_noback'
    
    PRESET_031 = '031_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target1on_noback'
    PRESET_032 = '032_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target2on_noback'
    PRESET_033 = '033_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target3on_noback'
    PRESET_034 = '034_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target4on_noback'
    PRESET_035 = '035_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target5on_noback'
    PRESET_036 = '036_geom_cell_cellstep_sheet_kcel_dock0456_kuhn088_target6on_noback'
    
    PRESET_037 = '037_geom_null_cellstep_linear_kcel_nodock_kuhn352_target1_noback'
    PRESET_038 = '038_geom_null_cellstep_linear_kcel_nodock_kuhn352_target2_noback'
    PRESET_039 = '039_geom_null_cellstep_linear_kcel_nodock_kuhn352_target3_noback'
    PRESET_040 = '040_geom_null_cellstep_linear_kcel_nodock_kuhn352_target4_noback'
    PRESET_041 = '041_geom_null_cellstep_linear_kcel_nodock_kuhn352_target5_noback'
    PRESET_042 = '042_geom_null_cellstep_linear_kcel_nodock_kuhn352_target6_noback'
    
    PRESET_043 = '043_geom_null_cellstep_linear_kcel_nodock_kuhn704_target1_noback'
    PRESET_044 = '044_geom_null_cellstep_linear_kcel_nodock_kuhn704_target2_noback'
    PRESET_045 = '045_geom_null_cellstep_linear_kcel_nodock_kuhn704_target3_noback'
    PRESET_046 = '046_geom_null_cellstep_linear_kcel_nodock_kuhn704_target4_noback'
    PRESET_047 = '047_geom_null_cellstep_linear_kcel_nodock_kuhn704_target5_noback'
    PRESET_048 = '048_geom_null_cellstep_linear_kcel_nodock_kuhn704_target6_noback'
    
    PRESET_049 = '049_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target1_noback'
    PRESET_050 = '050_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target2_noback'
    PRESET_051 = '051_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target3_noback'
    PRESET_052 = '052_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target4_noback'
    PRESET_053 = '053_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target5_noback'
    PRESET_054 = '054_geom_cell_cellstep_linear_kcel_nodock_kuhn176_target6_noback'
    
    PRESET_055 = '055_geom_cell_cellstep_linear_kcel_dock0246p80_kuhn088_target1on_noback'
    PRESET_056 = '056_geom_null_cellstep_sheet_kcel_nodock_kuhn088_targetAll_noback'
    PRESET_057 = '057_geom_null_myostep_linear_kcel_nodock_kuhn088_target1on_noback'
    
    PRESET_058 = '058_geom_null_cellstep_linear_kcel_nodock_kuhn088_target1on_noback'
    PRESET_059 = '059_geom_null_cell3step_linear_kcel_nodock_kuhn088_target1on_noback'
    PRESET_060 = '060_geom_null_cell5step_linear_kcel_nodock_kuhn088_target1on_noback'
    PRESET_061 = '061_geom_null_cell7step_linear_kcel_nodock_kuhn088_target1on_noback'
    
    def makeINAME(self,preset):
        self.IMAGE_NAME = preset + '.svg'
        self.IMAGE_JPG = preset + '.png'
        
        self.IMAGE_FULL_NAME = self.IMAGE_PATH + self.IMAGE_NAME
        self.IMAGE_FULL_JPG  = self.IMAGE_PATH + 'png/' + self.IMAGE_JPG
        
    def calcResiduesFromSegments(self,NS):
        # AA * CONTOUR / KUHN = NS
        # AA = NS * KUHN / CONTOUR
        
        if self.SET_AA is not None:
            return self.SET_AA
        
        return int(math.floor( NS * self.KUHN / self.CONTOUR ))
    
    def calcMaxSegments(self,AA = None):
        if (AA is None) :
            AA = self.DEFAULT_MAX_AA
        ms = int(math.floor( AA * self.CONTOUR / self.KUHN ))
        self.MAX_SEGMENT_NUM  = ms
        
        return ms
    
    def setPRESETGEOM(self,preset=None):
        if preset is None:
            preset = self.DEFAULT_PRESET_GEOM
        self.PRESETGEOM = preset
        
    def addKCAT(self):
        self.Kcat_TIMES = [
             0.53097345132743362831858407079646
            ,0.5
            ,0.38461538461538461538461538461538
            ,21.428571428571428571428571428571
            ,0.55555555555555555555555555555556
            ,0.35714285714285714285714285714286
            ,0.26905829596412556053811659192825
            ,16.216216216216216216216216216216
            ,5.8823529411764705882352941176471            
            ]
        
        self.Kcat_TIMES_LOWER = [
             0.26905829596412556053811659192825
            ,5.8823529411764705882352941176471
            ]
    
    def addLowerLL(self):
        self.Lower_LL = [35, 57]
        
    def addHigherLL(self):
        self.Higher_LL = [ ]
        
    def removeBackround(self):
        self.Kcat_TIMES = []
        self.Lower_LL   = []
        self.Higher_LL  = []
        self.Kcat_TIMES_LOWER = []
        
    def addBackground(self):
        self.addKCAT()
        self.addLowerLL()
        self.addHigherLL()
        
    def setSegment(self,segment=None,aa=None):
        if aa is not None :
            segment = math.ceil( float(aa) * self.CONTOUR / self.KUHN )
            self.SET_AA = aa
           
        self.MIN_SEGMENT_NUM = segment
        self.MAX_SEGMENT_NUM = segment + 1
            
        
    def doPRESET_001(self):
        self.makeINAME(self.PRESET_001)
        self.addBackground()
            
        self.KUHN = self.DEFAULT_KUHN
        
        self.MAX_SEGMENT_NUM = self.calcMaxSegments()
        if (self.USE_AA):
            self.XLIM_MAX = self.MAX_AA
        else :
            self.XLIM_MAX = self.MAX_SEGMENT_NUM
        self.STEP_X  = self.CELLULOSE_STEP_X
        self.STEP_Y  = self.CELLULOSE_STEP_Y
        
        self.MIN_TARGET_STEP_X = 1
        self.MAX_TARGET_STEP_X = 1
        self.MIN_BACK_STEP_X   = 0
        self.MAX_BACK_STEP_X   = None
        
        self.MIN_TARGET_STEP_Y = 0
        self.MAX_TARGET_STEP_Y = 0
        self.MIN_BACK_STEP_Y   = 0
        self.MAX_BACK_STEP_Y   = None
        
        self.KON   = self.CEL7A_KON
        self.KON_M = self.CEL7A_KON_M
        
        self.MIN_DOCKED     = 0
        self.MAX_DOCKED     = 10
        self.DOCK_STEP      = 4
        self.DOCK_STEP_NUM  = None   
        
        self.DOCKLIST = self.DOCKLIST_NODOCK
        
        self.PRESETGEOM = self.PRESETGEOM_NULL
        
        self.IS_SHEET_SUBSTRATE = False
        
        self.YLIM_MIN = 1e-5
        self.YLIM_MAX = 1e-1
        
        self.DRAW_HLINE = True
    
    def doPRESET_007(self):
        self.doPRESET_001()
        self.KUHN = 1.76
        self.MIN_SEGMENT_NUM = 4
        self.MAX_SEGMENT_NUM = self.calcMaxSegments()
        
        if (self.USE_AA):
            self.XLIM_MAX = self.MAX_AA
        else :
            self.XLIM_MAX = self.MAX_SEGMENT_NUM
    
    def doPRESET_013(self):
        self.doPRESET_001()
        self.MAX_TARGET_STEP_X = 1000
        
    def doPRESET_019(self):
        self.doPRESET_013()
        self.PRESETGEOM = self.PRESETGEOM_CELLULASE
        self.MAX_TARGET_STEP_X = self.MIN_TARGET_STEP_X
        
    def doPRESET_025(self):
        self.doPRESET_019()
        self.DOCKLIST = self.DOCKLIST_HIGH2
    
    def doPRESET_031(self):
        self.doPRESET_025()
        self.MAX_TARGET_STEP_X = 1000
        self.MIN_TARGET_STEP_X = -1000
        self.MIN_TARGET_STEP_Y = -1000
        self.MAX_TARGET_STEP_Y = 1000
        self.IS_SHEET_SUBSTRATE = True
    
    def doPRESET_037(self):
        self.doPRESET_001()
        self.KUHN = 3.52
        self.MIN_SEGMENT_NUM = 4
        self.MAX_SEGMENT_NUM = self.calcMaxSegments()
        
        if (self.USE_AA):
            self.XLIM_MAX = self.MAX_AA
        else :
            self.XLIM_MAX = self.MAX_SEGMENT_NUM
    
    def doPRESET_043(self):
        self.doPRESET_001()
        self.KUHN = 7.04
        self.MIN_SEGMENT_NUM = 4
        self.MAX_SEGMENT_NUM = self.calcMaxSegments()
        
        if (self.USE_AA):
            self.XLIM_MAX = self.MAX_AA
        else :
            self.XLIM_MAX = self.MAX_SEGMENT_NUM
    
    def doPRESET_049(self):
        self.doPRESET_007()
        self.PRESETGEOM = self.PRESETGEOM_CELLULASE
    
    def doPRESET_055(self):
        self.doPRESET_025()
        self.DOCKPERCENT = self.DOCKPERCENT_MID
        self.YLIM_MIN = 1e-3
        self.YLIM_MAX = 1e+1
    
    def doPRESET_056(self):
        self.doPRESET_001()
        self.MAX_TARGET_STEP_X = 1000
        self.MIN_TARGET_STEP_X = -1000
        self.MIN_TARGET_STEP_Y = -1000
        self.MAX_TARGET_STEP_Y = 1000
        self.IS_SHEET_SUBSTRATE = True
        
        self.YLIM_MIN = 5e-7
        self.YLIM_MAX = 1e-1
    
    def doPRESET_057(self):
        self.doPRESET_001()
        self.STEP_X  = self.MYOSIN_STEP_X
        self.STEP_Y  = self.MYOSIN_STEP_Y
        self.MAX_TARGET_STEP_X = 1000
        
    def applyPreset(self,preset):
        if preset == self.PRESET_001 :
            self.doPRESET_001()
            self.makeINAME(preset)
        elif preset == self.PRESET_002 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_003 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_004 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_005 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_006 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
        
        elif preset == self.PRESET_007 :
            self.doPRESET_007()
            self.makeINAME(preset)
        elif preset == self.PRESET_008 :
            self.doPRESET_007()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_009 :
            self.doPRESET_007()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_010 :
            self.doPRESET_007()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_011 :
            self.doPRESET_007()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_012 :
            self.doPRESET_007()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_013 :
            self.doPRESET_013()
            self.makeINAME(preset)
        elif preset == self.PRESET_014 :
            self.doPRESET_013()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
        elif preset == self.PRESET_015 :
            self.doPRESET_013()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
        elif preset == self.PRESET_016 :
            self.doPRESET_013()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
        elif preset == self.PRESET_017 :
            self.doPRESET_013()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
        elif preset == self.PRESET_018 :
            self.doPRESET_013()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_019 :
            self.doPRESET_019()
            self.makeINAME(preset)
        elif preset == self.PRESET_020 :
            self.doPRESET_019()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_021 :
            self.doPRESET_019()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_022 :
            self.doPRESET_019()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_023 :
            self.doPRESET_019()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_024 :
            self.doPRESET_019()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
        
        elif preset == self.PRESET_025 :
            self.doPRESET_025()
            self.makeINAME(preset)
        elif preset == self.PRESET_026 :
            self.doPRESET_025()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_027 :
            self.doPRESET_025()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_028 :
            self.doPRESET_025()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_029 :
            self.doPRESET_025()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_030 :
            self.doPRESET_025()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_031 :
            self.doPRESET_031()
            self.makeINAME(preset)
        elif preset == self.PRESET_032 :
            self.doPRESET_031()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_033 :
            self.doPRESET_031()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
            
#             self.MIN_SEGMENT_NUM = 10
#             self.MAX_SEGMENT_NUM = 10
            
        elif preset == self.PRESET_034 :
            self.doPRESET_031()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_035 :
            self.doPRESET_031()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_036 :
            self.doPRESET_031()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_037 :
            self.doPRESET_037()
            self.makeINAME(preset)
        elif preset == self.PRESET_038 :
            self.doPRESET_037()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_039 :
            self.doPRESET_037()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_040 :
            self.doPRESET_037()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_041 :
            self.doPRESET_037()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_042 :
            self.doPRESET_037()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_043 :
            self.doPRESET_043()
            self.makeINAME(preset)
        elif preset == self.PRESET_044 :
            self.doPRESET_043()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_045 :
            self.doPRESET_043()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_046 :
            self.doPRESET_043()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_047 :
            self.doPRESET_043()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_048 :
            self.doPRESET_043()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
            
        elif preset == self.PRESET_049 :
            self.doPRESET_049()
            self.makeINAME(preset)
        elif preset == self.PRESET_050 :
            self.doPRESET_049()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 2
            self.MAX_TARGET_STEP_X = 2
        elif preset == self.PRESET_051 :
            self.doPRESET_049()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 3
            self.MAX_TARGET_STEP_X = 3
        elif preset == self.PRESET_052 :
            self.doPRESET_049()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 4
            self.MAX_TARGET_STEP_X = 4
        elif preset == self.PRESET_053 :
            self.doPRESET_049()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 5
            self.MAX_TARGET_STEP_X = 5
        elif preset == self.PRESET_054 :
            self.doPRESET_049()
            self.makeINAME(preset)
            self.MIN_TARGET_STEP_X = 6
            self.MAX_TARGET_STEP_X = 6
        elif preset == self.PRESET_055 :
            self.doPRESET_055()
            self.makeINAME(preset)
        elif preset == self.PRESET_056 :
            self.doPRESET_056()
            self.makeINAME(preset)
        elif preset == self.PRESET_057 :
            self.doPRESET_057()
            self.makeINAME(preset)
            self.YLIM_MIN=None
            self.YLIM_MAX=None
        
        elif preset == self.PRESET_058 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.STEP_X  = self.CELLULOSE_STEP_X
            self.STEP_Y  = self.CELLULOSE_STEP_Y
            self.MIN_TARGET_STEP_X = 1
            self.MAX_TARGET_STEP_X = 1000
        
        elif preset == self.PRESET_059 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.STEP_X  = self.CELL_STEP3_X
            self.STEP_Y  = self.CELLULOSE_STEP_Y
            self.MIN_TARGET_STEP_X = 1
            self.MAX_TARGET_STEP_X = 1000
        
        elif preset == self.PRESET_060 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.STEP_X  = self.CELL_STEP5_X
            self.STEP_Y  = self.CELLULOSE_STEP_Y
            self.MIN_TARGET_STEP_X = 1
            self.MAX_TARGET_STEP_X = 1000
        
        elif preset == self.PRESET_061 :
            self.doPRESET_001()
            self.makeINAME(preset)
            self.STEP_X  = self.CELL_STEP7_X
            self.STEP_Y  = self.CELLULOSE_STEP_Y
            self.YLIM_MIN=3e-4
            self.YLIM_MAX=1e+0
            self.MIN_TARGET_STEP_X = 1
            self.MAX_TARGET_STEP_X = 1000
            
        
            
    
    def resetVars(self, preset = None):
        
        self.etol_abs = 1e-1
        self.etol_rel = 1e-1
        self.limit_num = 1
        self.cycle_num = 3
        self.option_dict = {'epsabs' : self.etol_abs, 'epsrel' : self.etol_rel, 'limit' : self.limit_num, 'limlst' : self.cycle_num}
        
        self.stepXNum = 80
        self.stepXHNum = 400
        self.stepYZNum = 6
        self.stepYZHNum = 24
        
        self.DOCKPERCENT = self.DEFAULT_DOCKPERCENT
        
        self.MIN_SEGMENT_NUM = 6
        self.SEGMENT_STEP = self.DEFAULT_SEGMENT_STEP
        
        self.MARKSIZE = 5
        
        self.YLIM_MIN = 1e-5
        self.YLIM_MAX = 1e-1
        
        self.XLIM_MIN = 0
        
        self.SET_AA = None
        
        self.removeBackround()
        
        if preset is None :
            preset = self.PRESET_001
    
        self.applyPreset(preset)
        
        
    def initVars(self, preset = None):
        self.resetVars(preset)
            
        if len(self.DOCKLIST) == 0 :
            self.DOCKLIST = list(range(self.MIN_DOCKED, self.MAX_DOCKED, self.DOCK_STEP))
        
        self.MAX_SEGMENT_NUM += 1
        
        self.MAX_TARGET_STEP_X += 1
        self.MAX_TARGET_STEP_Y += 1
        
        if self.MAX_BACK_STEP_X is None :
            self.MAX_BACK_STEP_X = self.MIN_BACK_STEP_X
        else :
            self.MAX_BACK_STEP_X += 1
            
        if self.MAX_BACK_STEP_Y is None :
            self.MAX_BACK_STEP_Y = self.MIN_BACK_STEP_Y
        else :
            self.MAX_BACK_STEP_Y += 1
        
        if self.MIN_TARGET_STEP_X < self.MAX_BACK_STEP_X :
            self.MIN_TARGET_STEP_X = self.MAX_BACK_STEP_X
        
        self.IMAGE_FULL_NAME = self.IMAGE_PATH + self.IMAGE_NAME
        self.IMAGE_FULL_JPG  = self.IMAGE_PATH + 'png/' + self.IMAGE_JPG
        
    
    def __init__(self):
        self.initVars()
        
p = pniniVars()

def getVars():
    return p