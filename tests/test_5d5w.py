import pynautilus as n

def main():
    # cnautilus -mtzin test_data/5d5w/hklout.mtz -seqin test_data/5d5w/5d5w.fasta
    # -pdbin test_data/5d5w/xyzout.pdb -colin-fo FP,SIGFP -colin-free FREE
    # -colin-fc FWT,PHWT -pdbout test_output/5d5w/pdbout.pdb
    # const std::string& seqin,
    #         const std::string& pdbin,
    #         const std::string& predin,
    #         const std::string& colin_fo,
    #         const std::string& colin_hl,
    #         const std::string& colin_phifom,
    #         const std::string& colin_fc,
    #         const std::string& colin_free

    input = n.Input(
        "test_data/5d5w/hklout.mtz",
        "test_data/5d5w/5d5w.fasta",
        "test_data/5d5w/xyzout.pdb",
        "test_data/5d5w/phosphate.map",
        "FP,SIGFP",
        "",
        "",
        "FWT,PHWT",
        "FREE"
    )

    output = n.Output("test_output/5d5w/pynautilus/pynautout.pdb")

    n.run(input, output, 3)

if __name__ == "__main__":
    main()