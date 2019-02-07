#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

void CPropertiesWnd::add_calc_ctrls()
{
  CMFCPropertyGridProperty* pCalcGroup = new CMFCPropertyGridProperty(_T("Calculators"));

  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    EvaporatingParticle::CCalculator* pCalc = pCalcColl->at(i);
    int nType = pCalc->type();
    switch(nType)
    {
      case EvaporatingParticle::CCalculator::ctPlaneYZ:
      {
        EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)pCalc;
    
        CMFCPropertyGridProperty* pXPlaneCalcGroup = new CMFCPropertyGridProperty(pPlaneYZCalc->get_name());

        CCalcResponseProperty* pXPlanePos = new CCalcResponseProperty(this, _T("Cross-Section X, mm"), COleVariant(10 * pPlaneYZCalc->get_plane_x()), _T("Set x-coordinate of the vertical cross-section plane."), pPlaneYZCalc->get_plane_x_ptr());
        pXPlaneCalcGroup->AddSubItem(pXPlanePos);

        CPlaneYZCalcCheckBox* pEnableCrossSectDrawing = new CPlaneYZCalcCheckBox(this, _T("Enable Cross-Section Drawing"), (_variant_t)pPlaneYZCalc->get_enable(), _T("Turns ON/OFF the cross-section drawing."), pPlaneYZCalc->get_enable_ptr());
        pXPlaneCalcGroup->AddSubItem(pEnableCrossSectDrawing);

        COleVariant var(pPlaneYZCalc->get_var_name(pPlaneYZCalc->get_clc_var_type()));
        CCalcResponseProperty* pVarSelectProp = new CCalcResponseProperty(this, _T("Variable to be Calculated"), var, _T("Specify a variable to be calculated at the selected cross-section."), pPlaneYZCalc->get_clc_var_type_ptr());
        for(int j = 0; j < EvaporatingParticle::CCalculator::clcCount; j++)
          pVarSelectProp->AddOption(pPlaneYZCalc->get_var_name(j));

        pVarSelectProp->AllowEdit(FALSE);
        pXPlaneCalcGroup->AddSubItem(pVarSelectProp);

// Character Length specially for Reynolds number calculations:
        if(pCalc->get_clc_var_type() == EvaporatingParticle::CCalculator::clcAveRe)
        {
          CMFCPropertyGridProperty* pCharLen = new CMFCPropertyGridProperty(_T("Character Length, mm"), COleVariant(10 * pCalc->get_char_length()), _T("Define the character length for Reynolds number calculation."), pCalc->get_char_length_ptr());
          pXPlaneCalcGroup->AddSubItem(pCharLen);
        }

        CString sPropTitle = pPlaneYZCalc->units();
        CMFCPropertyGridProperty* pCalcRes = new CMFCPropertyGridProperty(sPropTitle, COleVariant(pPlaneYZCalc->get_result()), _T("This read-only edit-line reads the result of calculation."), pPlaneYZCalc->get_result_ptr());
        pCalcRes->AllowEdit(FALSE);
        pXPlaneCalcGroup->AddSubItem(pCalcRes);

// Sequential calculations:
        CMFCPropertyGridProperty* pSeqCalcGroup = new CMFCPropertyGridProperty(_T("Sequential Calculations"));

        CMFCPropertyGridProperty* pStartPos = new CMFCPropertyGridProperty(_T("Start Position, mm"), COleVariant(10 * pPlaneYZCalc->get_start_pos()), _T("Set x-coordinate of the start position of the cross-section plane."), pPlaneYZCalc->get_start_pos_ptr());
        pSeqCalcGroup->AddSubItem(pStartPos);

        CMFCPropertyGridProperty* pStepX = new CMFCPropertyGridProperty(_T("End Position, mm"), COleVariant(10 * pPlaneYZCalc->get_end_pos()), _T("The position of the cross-section plane will be incremented by this step."), pPlaneYZCalc->get_end_pos_ptr());
        pSeqCalcGroup->AddSubItem(pStepX);

        CMFCPropertyGridProperty* pStepCount = new CMFCPropertyGridProperty(_T("Count of Steps"), COleVariant((long)pPlaneYZCalc->get_seq_calc_count()), _T("Set the total count of cross-section plane positions for this sequence."), pPlaneYZCalc->get_seq_calc_count_ptr());
        pSeqCalcGroup->AddSubItem(pStepCount);

// Start sequential calculations button:
        COleVariant var2(_T(""));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), var2, _T("Click to start a sequence of calculations."), (DWORD_PTR)pPlaneYZCalc);
        pSeqCalcGroup->AddSubItem(pStartCalcBtn);

// Output filename for sequential calculations:
        CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Output File"));
        static TCHAR BASED_CODE szFilter[] = _T("Data Files(*.csv)|*.csv|All Files(*.*)|*.*||");
        CMFCPropertyGridFileProperty* pOutFileProp = new CMFCPropertyGridFileProperty(_T("Select Filename"), TRUE, pPlaneYZCalc->get_filename(), _T("csv"),
          OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of the sequential calculations data file."), pPlaneYZCalc->get_filename_ptr());

        pDataGroup->AddSubItem(pOutFileProp);
        pSeqCalcGroup->AddSubItem(pDataGroup);

        pXPlaneCalcGroup->AddSubItem(pSeqCalcGroup);

// Remove calculator button:
        COleVariant var1(_T(""));
        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), var1, _T("Click to delete this calculator."), (DWORD_PTR)pPlaneYZCalc);
        pXPlaneCalcGroup->AddSubItem(pRemoveCalcBtn);

        pCalcGroup->AddSubItem(pXPlaneCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctSelRegions:
      {
        EvaporatingParticle::CSelectedRegionCalculator* pSelRegCalc = (EvaporatingParticle::CSelectedRegionCalculator*)pCalc;
    
        CMFCPropertyGridProperty* pSelRegCalcGroup = new CMFCPropertyGridProperty(pSelRegCalc->get_name());

// 2D regions selector:
        CString cRegNames = pSelRegCalc->get_sel_reg_names();
        CSelectRegionButton* pSelectRegButton = new CSelectRegionButton(this, _T("2D Regions"), cRegNames, _T("Click to select 2D regions at which the calculations will be performed."), pSelRegCalc->get_sel_reg_names_ptr());
        pSelRegCalcGroup->AddSubItem(pSelectRegButton);

        COleVariant var(pSelRegCalc->get_var_name(pSelRegCalc->get_clc_var_type()));
        CCalcResponseProperty* pVarSelectProp = new CCalcResponseProperty(this, _T("Variable to be Calculated"), var, _T("Specify a variable to be calculated at the selected cross-section."), pSelRegCalc->get_clc_var_type_ptr());
        for(int j = 0; j < EvaporatingParticle::CCalculator::clcCount; j++)
          pVarSelectProp->AddOption(pSelRegCalc->get_var_name(j));

        pVarSelectProp->AllowEdit(FALSE);
        pSelRegCalcGroup->AddSubItem(pVarSelectProp);

        CString sPropTitle = pSelRegCalc->units();
        CMFCPropertyGridProperty* pCalcRes = new CMFCPropertyGridProperty(sPropTitle, COleVariant(pSelRegCalc->get_result()), _T("This read-only edit-line reads the result of calculation."), pSelRegCalc->get_result_ptr());
        pCalcRes->AllowEdit(FALSE);
        pSelRegCalcGroup->AddSubItem(pCalcRes);

// Remove calculator button:
        COleVariant var1(_T(""));
        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), var1, _T("Click to delete this calculator."), (DWORD_PTR)pSelRegCalc);
        pSelRegCalcGroup->AddSubItem(pRemoveCalcBtn);

        pCalcGroup->AddSubItem(pSelRegCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongLine:
      {
        EvaporatingParticle::CLineCalculator* pLineCalc = (EvaporatingParticle::CLineCalculator*)pCalc;

        CMFCPropertyGridProperty* pLineCalcGroup = new CMFCPropertyGridProperty(pLineCalc->get_name());

        CMFCPropertyGridProperty* pProp;
        CMFCPropertyGridProperty* pStartPosGroup = new CMFCPropertyGridProperty(_T("Line Start Position, mm"), pLineCalc->get_start_ptr());
        pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pLineCalc->get_start().x), _T("X-coordinate of the line starting point."));
        pStartPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pLineCalc->get_start().y), _T("Y-coordinate of the line starting point."));
        pStartPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pLineCalc->get_start().z), _T("Z-coordinate of the line starting point."));
        pStartPosGroup->AddSubItem(pProp);
        pLineCalcGroup->AddSubItem(pStartPosGroup);

        CMFCPropertyGridProperty* pEndPosGroup = new CMFCPropertyGridProperty(_T("Line End Position, mm"), pLineCalc->get_end_ptr());
        pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pLineCalc->get_end().x), _T("X-coordinate of the line starting point."));
        pEndPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pLineCalc->get_end().y), _T("Y-coordinate of the line starting point."));
        pEndPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pLineCalc->get_end().z), _T("Z-coordinate of the line starting point."));
        pEndPosGroup->AddSubItem(pProp);
        pLineCalcGroup->AddSubItem(pEndPosGroup);

        CMFCPropertyGridProperty* pStepCount = new CMFCPropertyGridProperty(_T("Count of Steps"), COleVariant((long)pLineCalc->get_steps_count()), _T("Set the count of steps along the line."), pLineCalc->get_steps_count_ptr());
        pLineCalcGroup->AddSubItem(pStepCount);

        COleVariant var(pLineCalc->get_var_name(pLineCalc->get_clc_var_type()));
        CCalcResponseProperty* pVarSelectProp = new CCalcResponseProperty(this, _T("Variable to be Calculated"), var, _T("Specify a variable to be calculated at the selected cross-section."), pLineCalc->get_clc_var_type_ptr());
        for(int j = 0; j < pLineCalc->calc_vars_count(); j++)
          pVarSelectProp->AddOption(pLineCalc->get_var_name(j));

        pVarSelectProp->AllowEdit(FALSE);
        pLineCalcGroup->AddSubItem(pVarSelectProp);

// Character Length specially for Reynolds number calculations:
        if(pCalc->get_clc_var_type() == EvaporatingParticle::CLineCalculator::lcRe)
        {
          CMFCPropertyGridProperty* pCharLen = new CMFCPropertyGridProperty(_T("Character Length, mm"), COleVariant(10 * pCalc->get_char_length()), _T("Define the character length for Reynolds number calculation."), pCalc->get_char_length_ptr());
          pLineCalcGroup->AddSubItem(pCharLen);
        }

// Start calculations button:
        COleVariant var2(_T(""));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), var2, _T("Click to start calculations."), (DWORD_PTR)pLineCalc);
        pLineCalcGroup->AddSubItem(pStartCalcBtn);

// Output filename for line calculator:
        CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Output File"));
        static TCHAR BASED_CODE szFilter[] = _T("Data Files(*.csv)|*.csv|All Files(*.*)|*.*||");
        CMFCPropertyGridFileProperty* pOutFileProp = new CMFCPropertyGridFileProperty(_T("Select Filename"), TRUE, pLineCalc->get_filename(), _T("csv"),
          OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location for the result data file."), pLineCalc->get_filename_ptr());

        pDataGroup->AddSubItem(pOutFileProp);
        pLineCalcGroup->AddSubItem(pDataGroup);

// Remove calculator button:
        COleVariant var1(_T(""));
        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), var1, _T("Click to delete this calculator."), (DWORD_PTR)pLineCalc);
        pLineCalcGroup->AddSubItem(pRemoveCalcBtn);

        pCalcGroup->AddSubItem(pLineCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCalc:
      {
        EvaporatingParticle::CTrackCalculator* pTrackCalc = (EvaporatingParticle::CTrackCalculator*)pCalc;

        CMFCPropertyGridProperty* pTrackCalcGroup = new CMFCPropertyGridProperty(pTrackCalc->get_name());

        CCalcResponseProperty* pXPlanePos = new CCalcResponseProperty(this, _T("Cross-Section X, mm"), COleVariant(10 * pTrackCalc->get_cs_pos()), _T("Set x-coordinate of the vertical cross-section plane."), pTrackCalc->get_cs_pos_ptr());
        pTrackCalcGroup->AddSubItem(pXPlanePos);

        COleVariant var(pTrackCalc->get_var_name(pTrackCalc->get_clc_var_type()));
        CCalcResponseProperty* pVarSelectProp = new CCalcResponseProperty(this, _T("Variable to be Calculated"), var, _T("Specify a variable to be calculated at the selected cross-section."), pTrackCalc->get_clc_var_type_ptr());
        for(int j = 0; j < EvaporatingParticle::CTrackCalculator::clcTrackCount; j++)
          pVarSelectProp->AddOption(pTrackCalc->get_var_name(j));

        pVarSelectProp->AllowEdit(FALSE);
        pTrackCalcGroup->AddSubItem(pVarSelectProp);

        CString sPropTitle = pTrackCalc->units();
        CMFCPropertyGridProperty* pCalcRes = new CMFCPropertyGridProperty(sPropTitle, COleVariant(pTrackCalc->get_result()), _T("This read-only edit-line reads the result of calculation."), pTrackCalc->get_result_ptr());
        pCalcRes->AllowEdit(FALSE);
        pTrackCalcGroup->AddSubItem(pCalcRes);

// Sequential calculations:
        CMFCPropertyGridProperty* pSeqCalcGroup = new CMFCPropertyGridProperty(_T("Sequential Calculations"));

        CMFCPropertyGridProperty* pStartPos = new CMFCPropertyGridProperty(_T("Start X, mm"), COleVariant(10 * pTrackCalc->get_start_x()), _T("Set x-coordinate of the start cross-section."), pTrackCalc->get_start_x_ptr());
        pSeqCalcGroup->AddSubItem(pStartPos);

        CMFCPropertyGridProperty* pEndPos = new CMFCPropertyGridProperty(_T("End X, mm"), COleVariant(10 * pTrackCalc->get_end_x()), _T("Set x-coordinate of the end cross-section."), pTrackCalc->get_end_x_ptr());
        pSeqCalcGroup->AddSubItem(pEndPos);

        CMFCPropertyGridProperty* pCSCount = new CMFCPropertyGridProperty(_T("Count of Cross-Sections"), COleVariant((long)pTrackCalc->get_cs_count()), _T("Set the total count of cross-sections."), pTrackCalc->get_cs_count_ptr());
        pSeqCalcGroup->AddSubItem(pCSCount);

// Start calculations button:
        COleVariant var3(_T(""));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), var3, _T("Click to start calculations."), (DWORD_PTR)pTrackCalc);
        pSeqCalcGroup->AddSubItem(pStartCalcBtn);

// Output filename for the track calculator:
        CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Output File"));
        static TCHAR BASED_CODE szFilter[] = _T("Data Files(*.csv)|*.csv|All Files(*.*)|*.*||");
        CMFCPropertyGridFileProperty* pOutFileProp = new CMFCPropertyGridFileProperty(_T("Select Filename"), TRUE, pTrackCalc->get_filename(), _T("csv"),
          OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location for the result data file."), pTrackCalc->get_filename_ptr());

        pDataGroup->AddSubItem(pOutFileProp);
        pSeqCalcGroup->AddSubItem(pDataGroup);
        pTrackCalcGroup->AddSubItem(pSeqCalcGroup);
// Remove calculator button:
        COleVariant var1(_T(""));
        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), var1, _T("Click to delete this calculator."), (DWORD_PTR)pTrackCalc);
        pTrackCalcGroup->AddSubItem(pRemoveCalcBtn);

        pCalcGroup->AddSubItem(pTrackCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCrossSect:
      {
        EvaporatingParticle::CTrackCrossSectionCalculator* pTrackCalc = (EvaporatingParticle::CTrackCrossSectionCalculator*)pCalc;
        EvaporatingParticle::CColoredCrossSection* pObj = pTrackCalc->get_object();

        CMFCPropertyGridProperty* pTrackCalcGroup = new CMFCPropertyGridProperty(pTrackCalc->get_name());

        CMFCPropertyGridProperty* pXPlanePos = new CMFCPropertyGridProperty(_T("Cross-Section X, mm"), COleVariant(10 * pObj->get_cross_sect_pos()), _T("Set x-coordinate of the vertical cross-section plane."), pObj->get_cross_sect_pos_ptr());
        pTrackCalcGroup->AddSubItem(pXPlanePos);

        COleVariant var(pObj->get_var_name(pObj->get_var()));
        CMFCPropertyGridProperty* pVarSelectProp = new CMFCPropertyGridProperty(_T("Variable"), var, _T("Specify the variable, value of which will be used to color the points in the cross-section."), pObj->get_var_ptr());
        for(int j = 0; j < EvaporatingParticle::CColoredCrossSection::varCount; j++)
          pVarSelectProp->AddOption(pObj->get_var_name(j));

        pVarSelectProp->AllowEdit(FALSE);
        pTrackCalcGroup->AddSubItem(pVarSelectProp);

// Start calculations button:
        COleVariant var3(_T(""));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), var3, _T("Click to start calculations."), (DWORD_PTR)pTrackCalc);
        pTrackCalcGroup->AddSubItem(pStartCalcBtn);

// Statistics, informational edit-lines:
        CMFCPropertyGridProperty* pStatisticsGroup = new CMFCPropertyGridProperty("Statistics");

        CMFCPropertyGridProperty* pPointsCount = new CMFCPropertyGridProperty(_T("Points Count"), COleVariant((long)pTrackCalc->get_count()), _T("Info: Count of points in the cross-section."), NULL);
        pPointsCount->AllowEdit(FALSE);
        pStatisticsGroup->AddSubItem(pPointsCount);
        CMFCPropertyGridProperty* pAverageGroup = new CMFCPropertyGridProperty("Average Coordinates");
        pStatisticsGroup->AddSubItem(pAverageGroup);

        CMFCPropertyGridProperty* pAverY = new CMFCPropertyGridProperty(_T("Y, mm"), COleVariant(10 * pTrackCalc->get_center().y), _T("Info: Average Y coordinate in the cross-section."), NULL);
        pAverY->AllowEdit(FALSE);
        pAverageGroup->AddSubItem(pAverY);
        CMFCPropertyGridProperty* pAverZ = new CMFCPropertyGridProperty(_T("Z, mm"), COleVariant(10 * pTrackCalc->get_center().z), _T("Info: Average Z coordinate in the cross-section."), NULL);
        pAverZ->AllowEdit(FALSE);
        pAverageGroup->AddSubItem(pAverZ);

        CMFCPropertyGridProperty* pSigmaGroup = new CMFCPropertyGridProperty("Mean Square Deviation");
        pStatisticsGroup->AddSubItem(pSigmaGroup);

        CMFCPropertyGridProperty* pSigmaY = new CMFCPropertyGridProperty(_T("Sigma Y, mm"), COleVariant(10 * pTrackCalc->get_sigma().y), _T("Info: Mean square deviation in Y direction in the cross-section."), NULL);
        pSigmaY->AllowEdit(FALSE);
        pSigmaGroup->AddSubItem(pSigmaY);
        CMFCPropertyGridProperty* pSigmaZ = new CMFCPropertyGridProperty(_T("Sigma Z, mm"), COleVariant(10 * pTrackCalc->get_sigma().z), _T("Info: Mean square deviation in Z direction in the cross-section."), NULL);
        pSigmaZ->AllowEdit(FALSE);
        pSigmaGroup->AddSubItem(pSigmaZ);

        pTrackCalcGroup->AddSubItem(pStatisticsGroup);

// Output filename for the cross-section track calculator:
        CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Output File"));
        static TCHAR BASED_CODE szFilter[] = _T("Data Files(*.csv)|*.csv|All Files(*.*)|*.*||");
        CMFCPropertyGridFileProperty* pOutFileProp = new CMFCPropertyGridFileProperty(_T("Select Filename"), TRUE, pTrackCalc->get_filename(), _T("csv"),
          OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location for the result data file."), pTrackCalc->get_filename_ptr());

        pDataGroup->AddSubItem(pOutFileProp);
        pTrackCalcGroup->AddSubItem(pDataGroup);

// Remove calculator button:
        COleVariant var1(_T(""));
        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), var1, _T("Click to delete this calculator."), (DWORD_PTR)pTrackCalc);
        pTrackCalcGroup->AddSubItem(pRemoveCalcBtn);

        pCalcGroup->AddSubItem(pTrackCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongSelTracks:
      {
        EvaporatingParticle::CSelectedTracksCalculator* pSelTrackCalc = (EvaporatingParticle::CSelectedTracksCalculator*)pCalc;
        CMFCPropertyGridProperty* pTrackCalcGroup = new CMFCPropertyGridProperty(pSelTrackCalc->get_name());

        CMFCPropertyGridProperty* pPointsCount = new CMFCPropertyGridProperty(_T("Skipped Points Count"), COleVariant((long)pSelTrackCalc->get_skip_points_count()), _T("Every n-th point in the track will be skipped."), pSelTrackCalc->get_skip_points_count_ptr());
        pTrackCalcGroup->AddSubItem(pPointsCount);

        CSelectTrajectButton* pSelTrackBtn = new CSelectTrajectButton(this, _T("Trajectory Selection"), _T(""), _T("Click this button to enter/exit the trajectory selection tool."), NULL);
        pSelTrackBtn->SetValue(pSelTrackBtn->ButtonValue());
        pTrackCalcGroup->AddSubItem(pSelTrackBtn);

        COleVariant sDefFolder(pSelTrackCalc->get_out_folder());
        CSelectFolderButton* pOutFolderBtn = new CSelectFolderButton(this, _T("Output Folder"), sDefFolder, _T("Click and browse to select the output folder."), (DWORD_PTR)pSelTrackCalc);
        pTrackCalcGroup->AddSubItem(pOutFolderBtn);

// Enable/disable different forces calculation:
        CMFCPropertyGridProperty* pEnableGroup = new CMFCPropertyGridProperty(_T("Select forces to calculate"));
        CCheckBoxButton* pBtn = new CCheckBoxButton(this, _T("Gas Drag"), _T(""), _T("Enable/disable gas drag force calculation."), (DWORD_PTR)pSelTrackCalc->get_enable_gas_drag_ptr());
        pEnableGroup->AddSubItem(pBtn);

        pBtn = new CCheckBoxButton(this, _T("External DC Field"), _T(""), _T("Enable/disable external DC field calculation."), (DWORD_PTR)pSelTrackCalc->get_enable_dc_field_ptr());
        pEnableGroup->AddSubItem(pBtn);

        pBtn = new CCheckBoxButton(this, _T("External RF Field"), _T(""), _T("Enable/disable external RF field calculation."), (DWORD_PTR)pSelTrackCalc->get_enable_rf_field_ptr());
        pEnableGroup->AddSubItem(pBtn);

        pBtn = new CCheckBoxButton(this, _T("Space Charge Field"), _T(""), _T("Enable/disable space charge field calculation."), (DWORD_PTR)pSelTrackCalc->get_enable_clmb_ptr());
        pEnableGroup->AddSubItem(pBtn);

        pTrackCalcGroup->AddSubItem(pEnableGroup);

// Actions:
        EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
        EvaporatingParticle::CIdsVector vIds = pDrawObj->get_sel_traject_ids();
        bool bEnableStart = vIds.size() > 0;
        CMFCPropertyGridProperty* pActionsGroup = new CMFCPropertyGridProperty(_T("Actions"));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), _T(""), _T("Click to start calculations."), (DWORD_PTR)pSelTrackCalc);
        pStartCalcBtn->Enable(bEnableStart);
        pActionsGroup->AddSubItem(pStartCalcBtn);

        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), _T(""), _T("Click to delete this calculator."), (DWORD_PTR)pSelTrackCalc);
        pActionsGroup->AddSubItem(pRemoveCalcBtn);

        pTrackCalcGroup->AddSubItem(pActionsGroup);

        pCalcGroup->AddSubItem(pTrackCalcGroup);
        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackFaceCross:
      {
        EvaporatingParticle::CTackFaceCross* pCrossFacesCalc = (EvaporatingParticle::CTackFaceCross*)pCalc;
        CMFCPropertyGridProperty* pFaceCalcGroup = new CMFCPropertyGridProperty(pCrossFacesCalc->get_name());

        CMFCPropertyGridProperty* pStartPos = new CMFCPropertyGridProperty(_T("Start X, mm"), COleVariant(10 * pCrossFacesCalc->get_start_coord()), _T("Set the X-coordinate of the lower bound of the selection region."), pCrossFacesCalc->get_start_coord_ptr());
        pFaceCalcGroup->AddSubItem(pStartPos);

        CMFCPropertyGridProperty* pEndPos = new CMFCPropertyGridProperty(_T("End X, mm"), COleVariant(10 * pCrossFacesCalc->get_end_coord()), _T("Set the X-coordinate of the upper bound of the selection region."), pCrossFacesCalc->get_end_coord_ptr());
        pFaceCalcGroup->AddSubItem(pEndPos);

        CMFCPropertyGridProperty* pActionsGroup = new CMFCPropertyGridProperty(_T("Actions"));
        CStartCalcButton* pStartCalcBtn = new CStartCalcButton(this, _T("Start Calculations"), _T(""), _T("Click to start calculations."), (DWORD_PTR)pCrossFacesCalc);
        pActionsGroup->AddSubItem(pStartCalcBtn);

        CRemoveCalcButton* pRemoveCalcBtn = new CRemoveCalcButton(this, _T("Remove Calculator"), _T(""), _T("Click to delete this calculator."), (DWORD_PTR)pCrossFacesCalc);
        pActionsGroup->AddSubItem(pRemoveCalcBtn);

        pFaceCalcGroup->AddSubItem(pActionsGroup);

        pCalcGroup->AddSubItem(pFaceCalcGroup);
        break;
      }
    }
  }

  m_wndPropList.AddProperty(pCalcGroup);
}

void CPropertiesWnd::set_calc_data()
{
  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    EvaporatingParticle::CCalculator* pCalc = pCalcColl->at(i);
    int nType = pCalc->type();
    switch(nType)
    {
      case EvaporatingParticle::CCalculator::ctPlaneYZ: 
      {
        EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)pCalc;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_start_pos_ptr());
        if(pProp != NULL)
          pPlaneYZCalc->set_start_pos(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_end_pos_ptr());
        if(pProp != NULL)
          pPlaneYZCalc->set_end_pos(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_seq_calc_count_ptr());
        if(pProp != NULL)
          pPlaneYZCalc->set_seq_calc_count(pProp->GetValue().lVal);
        
        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_filename_ptr());
        if(pProp != NULL)
        {
          CString cFile = (CString)pProp->GetValue();
          std::string str = std::string(CT2CA(cFile));
          pPlaneYZCalc->set_filename(str.c_str());
        }

        pProp = m_wndPropList.FindItemByData(pCalc->get_char_length_ptr());
        if(pProp != NULL)
          pCalc->set_char_length(0.1 * pProp->GetValue().dblVal);

        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongLine: 
      {
        EvaporatingParticle::CLineCalculator* pLineCalc = (EvaporatingParticle::CLineCalculator*)pCalc;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pLineCalc->get_start_ptr());
        if(pProp != NULL)
        {
          EvaporatingParticle::Vector3D vStartPos;
          vStartPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
          vStartPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
          vStartPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
          pLineCalc->set_start(vStartPos);
        }

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_end_ptr());
        if(pProp != NULL)
        {
          EvaporatingParticle::Vector3D vEndPos;
          vEndPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
          vEndPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
          vEndPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
          pLineCalc->set_end(vEndPos);
        }

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_steps_count_ptr());
        if(pProp != NULL)
          pLineCalc->set_steps_count(pProp->GetValue().lVal);

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_filename_ptr());
        if(pProp != NULL)
        {
          CString cFile = (CString)pProp->GetValue();
          std::string str = std::string(CT2CA(cFile));
          pLineCalc->set_filename(str.c_str());
        }

        pProp = m_wndPropList.FindItemByData(pCalc->get_char_length_ptr());
        if(pProp != NULL)
          pCalc->set_char_length(0.1 * pProp->GetValue().dblVal);

        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCalc:
      {
        EvaporatingParticle::CTrackCalculator* pTrackCalc = (EvaporatingParticle::CTrackCalculator*)pCalc;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pTrackCalc->get_start_x_ptr());
        if(pProp != NULL)
          pTrackCalc->set_start_x(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_end_x_ptr());
        if(pProp != NULL)
          pTrackCalc->set_end_x(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_cs_count_ptr());
        if(pProp != NULL)
          pTrackCalc->set_cs_count(pProp->GetValue().lVal);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_filename_ptr());
        if(pProp != NULL)
        {
          CString cFile = (CString)pProp->GetValue();
          std::string str = std::string(CT2CA(cFile));
          pTrackCalc->set_filename(str.c_str());
        }

        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCrossSect:
      {
        EvaporatingParticle::CTrackCrossSectionCalculator* pTrackCalc = (EvaporatingParticle::CTrackCrossSectionCalculator*)pCalc;
        EvaporatingParticle::CColoredCrossSection* pObj = pTrackCalc->get_object();

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_cross_sect_pos_ptr());
        if(pProp != NULL)
          pObj->set_cross_sect_pos(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pObj->get_var_ptr());
        if(pProp != NULL)
        {
          CString cName = (CString)pProp->GetValue();
          for(int j = 0; j < EvaporatingParticle::CColoredCrossSection::varCount; j++)
          {
            if(cName == pObj->get_var_name(j))
            {
              pObj->set_var(j);
              break;
            }
          }
        }

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_filename_ptr());
        if(pProp != NULL)
        {
          CString cFile = (CString)pProp->GetValue();
          std::string str = std::string(CT2CA(cFile));
          pTrackCalc->set_filename(str.c_str());
        }

        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongSelTracks:
      {
        EvaporatingParticle::CSelectedTracksCalculator* pSelTrackCalc = (EvaporatingParticle::CSelectedTracksCalculator*)pCalc;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSelTrackCalc->get_skip_points_count_ptr());
        if(pProp != NULL)
        {
          pSelTrackCalc->set_skip_points_count((UINT)pProp->GetValue().llVal);
        }

        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackFaceCross:
      {
        EvaporatingParticle::CTackFaceCross* pCrossFacesCalc = (EvaporatingParticle::CTackFaceCross*)pCalc;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pCrossFacesCalc->get_start_coord_ptr());
        if(pProp != NULL)
          pCrossFacesCalc->set_start_coord(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pCrossFacesCalc->get_end_coord_ptr());
        if(pProp != NULL)
          pCrossFacesCalc->set_end_coord(0.1 * pProp->GetValue().dblVal);

        break;
      }
    }
  }
}

void CPropertiesWnd::update_calc_ctrls()
{
  bool bReady = !m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready();

  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    EvaporatingParticle::CCalculator* pCalc = pCalcColl->at(i);
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pCalc->get_enable_ptr());
    if(pProp != NULL)
      pProp->Enable(bReady);

    int nType = pCalc->type();
    switch(nType)
    {
      case EvaporatingParticle::CCalculator::ctPlaneYZ: 
      {
        EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)pCalc;

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_plane_x_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_clc_var_type_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_result_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_start_pos_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_end_pos_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_seq_calc_count_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);
        
        pProp = m_wndPropList.FindItemByData(pPlaneYZCalc->get_filename_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        break;
      }
      case EvaporatingParticle::CCalculator::ctSelRegions:
      {
        EvaporatingParticle::CSelectedRegionCalculator* pSelRegCalc = (EvaporatingParticle::CSelectedRegionCalculator*)pCalc;

        pProp = m_wndPropList.FindItemByData(pSelRegCalc->get_sel_reg_names_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongLine:
      {
        EvaporatingParticle::CLineCalculator* pLineCalc = (EvaporatingParticle::CLineCalculator*)pCalc;

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_start_ptr());
        if(pProp != NULL)
        {
          pProp->Enable(bReady);
          pProp->GetSubItem(0)->Enable(bReady);
          pProp->GetSubItem(1)->Enable(bReady);
          pProp->GetSubItem(2)->Enable(bReady);
        }

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_end_ptr());
        if(pProp != NULL)
        {
          pProp->Enable(bReady);
          pProp->GetSubItem(0)->Enable(bReady);
          pProp->GetSubItem(1)->Enable(bReady);
          pProp->GetSubItem(2)->Enable(bReady);
        }

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_steps_count_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_clc_var_type_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pLineCalc->get_filename_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCalc:
      {
        EvaporatingParticle::CTrackCalculator* pTrackCalc = (EvaporatingParticle::CTrackCalculator*)pCalc;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pTrackCalc->get_cs_pos_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_start_x_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_end_x_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_cs_count_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        pProp = m_wndPropList.FindItemByData(pTrackCalc->get_filename_ptr());
        if(pProp != NULL)
          pProp->Enable(bReady);

        break;
      }
    }
  }
}
