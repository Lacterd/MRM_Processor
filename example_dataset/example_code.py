"""

Example command-line workflow for MRM Processor.

This script reproduces the RT correction example

presented in the manuscript.

"""
from MRM_Processor_Analyze_GUI import MRMProcess

def main():
    # STD data
    BA_STD = MRMProcess('BA-STD.mzML')
    BA_STD.load_TargetList('MRM_List_ref.xlsx','IS')
    BA_STD.detect_Peak_MRM('Analyze')
    
    # STD RT Shift data without RT correction
    BA_STD_RT_Shift_without_correction = MRMProcess('BA-STD-RT-Shift.mzML')
    BA_STD_RT_Shift_without_correction.load_TargetList_without_correction(
        'MRM_List_ref.xlsx',
        'IS_RT_Shift_without_correction'
        )
    BA_STD_RT_Shift_without_correction.detect_Peak_MRM(
        'Analyze_without_correction'
        )
    
    # STD RT Shift data with RT correction
    BA_STD_RT_Shift_with_correction = MRMProcess('BA-STD-RT-Shift.mzML')
    BA_STD_RT_Shift_with_correction.load_TargetList(
        'MRM_List_ref.xlsx',
        'IS_RT_Shift_with_correction'
        )
    BA_STD_RT_Shift_with_correction.detect_Peak_MRM(
        'Analyze_with_correction'
        )
    
    print("Example analysis completed successfully.")

if __name__ == "__main__":

    main()