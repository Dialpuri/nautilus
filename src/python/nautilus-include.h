#include <string>
#include <optional>

class NautilusInput { 
public:
    NautilusInput(
        const std::string& mtzin,
        const std::string& seqin,
        const std::string& pdbin,
        const std::string& predin, 
        const std::string& colin_fo, 
        const std::string& colin_hl, 
        const std::string& colin_phifom,
        const std::string& colin_fc,
        const std::string& colin_free
    ) {
        this->mtzin = mtzin;
        this->seqin = seqin;
        this->pdbin = pdbin;
        this->predin = predin; 
        this->colin_fo = colin_fo; 
        this->colin_fc = colin_fc; 
        this->colin_phifom = colin_phifom;
        this->colin_fc = colin_fc; 
        this->colin_free = colin_free;

        if (mtzin == "") { throw std::runtime_error("MTZ Path must not be empty");}
        if (seqin == "") { throw std::runtime_error("SEQ Path must not be empty");}
        if (pdbin == "") { throw std::runtime_error("PDB Path must not be empty");}
        if (predin == "") { throw std::runtime_error("Predicition Path must not be empty");}

    };

    std::string get_mtz_path() const { return mtzin; }
    std::string get_seq_path() const { return seqin; }
    std::string get_pdb_path() const { return pdbin; }
    std::string get_prediction_path() const { return predin; }

    std::optional<std::string> get_fobs() const { 
        if (colin_fo == "") {
            return std::nullopt;
        }
        return colin_fo;
    }

    std::optional<std::string> get_hl() const { 
        if (colin_hl == "") {
            return std::nullopt;
        }
        return colin_hl;
    }

    std::optional<std::string> get_phifom() const { 
        if (colin_phifom == "") {
            return std::nullopt;
        }
        return colin_phifom;
    }

    std::optional<std::string> get_fc() const { 
        if (colin_fc == "") {
            return std::nullopt;
        }
        return colin_fc;
    }

    std::optional<std::string> get_free() const { 
        if (colin_free == "") {
            return std::nullopt;
        }
        return colin_free;
    }

private:
    std::string mtzin;
    std::string seqin;
    std::string pdbin;
    std::string predin; 
    std::string colin_fo; 
    std::string colin_hl; 
    std::string colin_phifom;
    std::string colin_fc;
    std::string colin_free;
};

class NautilusOutput { 
public: 
    NautilusOutput(const std::string& pdbout) { 
        this->pdbout = pdbout; 

        if (pdbout == "") {throw std::runtime_error("PDB Out must not be empty");}
    }; 

    std::string get_pdb_out() const { return pdbout; }
private:
    std::string pdbout; 
};