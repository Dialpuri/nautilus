import { useEffect, useRef, useState } from "react";
import nautilus_module from "../../../wasm/pnautilus.js"
import {MoorhenContainer, MoorhenContextProvider, MoorhenMap, MoorhenMolecule} from 'moorhen'

export default function Index() { 

    const [selectedFile, setSelectedFile] = useState(null);
    const [cootInitialized, setCootInitialized] = useState(false)
    const controls = useRef();

    const [findState, setFindState] = useState(false);
    const [growState, setGrowState] = useState(false);
    const [joinState, setJoinState] = useState(false);
    const [linkState, setLinkState] = useState(false);
    const [pruneState, setPruneState] = useState(false);
    const [rebuildChainState, setRebuildChainState] = useState(false);
    const [sequenceState, setSequenceState] = useState(false);
    const [rebuildBasesState, setRebuildBasesState] = useState(false);

    
    const display_map = (map_data) => { 
        if (cootInitialized) {
            const map = new MoorhenMap(controls.current.commandCentre, controls.current.glRef);
            const mapMetadata = { 
                F: "FWT", 
                PHI: "PHWT", 
                Fobs: "FP", 
                SigFobs: "SIGFP", 
                FreeR: "FREE", 
                isDifference: false,
                useWeight: false, 
                calcStructFact: true,
                }
            map.loadToCootFromMtzData(map_data, "map-name", mapMetadata);
            controls.current.changeMaps({ action: "Add", item: map })
            controls.current.setActiveMap(map)
        }
    }

    const display_molecule = (data, molecule, id) => { 
        if (cootInitialized) { 
            let newMolecule = new MoorhenMolecule(controls.current.commandCentre, controls.current.glRef, controls.current.monomerLibrary)
            newMolecule.loadToCootFromString(data, id).then(() => {
                controls.current.changeMolecules({action: 'Add', item: newMolecule});
                newMolecule.fetchIfDirtyAndDraw('CAs').then(() => {
                        newMolecule.centreOn()
                    }
                )
            })

            molecule = newMolecule
        }
    }

    
    const forwardControls = (forwardedControls) => {
        setCootInitialized(true)
        controls.current = forwardedControls
    }
    

    let reader = new FileReader();
    let fileName = ""
    let fobs = "FP,SIGFP"
    let fcalc = "FWT,PHWT"

    const load_file = (file) => { 
        fileName = file.name
        reader.addEventListener('loadend', read_file);
        reader.readAsArrayBuffer(file);
    }

    const read_file = () => {
        console.log(reader, fileName)
        const map_data = new Uint8Array(reader.result);
        display_map(map_data)

        nautilus_module().then((Module) => {
            Module['FS_createDataFile']('/', fileName, map_data, true, true, true)
            Module.load_mtz(fileName, fobs, fcalc);

            var find = Module.FS.readFile("find.pdb", {encoding:'utf8'})
            display_molecule(find, setFindMolecule, 'find-1')
            
            var grow = Module.FS.readFile("grow.pdb", {encoding:'utf8'})
            display_molecule(grow, setGrowMolecule, 'grow-1')

            var join = Module.FS.readFile("join.pdb", {encoding:'utf8'})
            display_molecule(join, setJoinMolecule, 'join-1')

            var link = Module.FS.readFile("link.pdb", {encoding:'utf8'})
            display_molecule(link, setLinkMolecule, 'link-1')

            var prune = Module.FS.readFile("prune.pdb", {encoding:'utf8'})
            display_molecule(prune, setPruneMolecule, 'prune-1')

            var rebuild_chain = Module.FS.readFile("rebuild_chain.pdb", {encoding:'utf8'})
            display_molecule(rebuild_chain, setRebuildChainMolecule, 'rebuild_chain-1')

            var sequence = Module.FS.readFile("sequence.pdb", {encoding:'utf8'})
            display_molecule(sequence, setSequenceMolecule, 'sequence-1')

            var rebuild_bases = Module.FS.readFile("rebuild_bases.pdb", {encoding:'utf8'})
            display_molecule(rebuild_bases, setRebuildBasesMolecule, 'rebuild_bases-1')

        })

    };

    return (
        <div className="flex flex-col align-middle justify-center">  

        <label for="colin-fo">F_obs column labels</label>      
        <input type="text" id="colin-fo" value="FP,SIGFP" onChange={(e) => {fobs=e.target.value}}/>

        <label for="colin-fc">F_calc column labels</label>      
        <input type="text" id="colin-fc" value="FWT,PHWT" onChange={(e) => {fcalc=e.target.value}}/>

        <input type="file" onChange={(e) => load_file(e.target.files[0])}/>

        {/* <div className="flex flex-wrap justify-center">
            <button onClick={() => {findState ? setFindState(false): setFindState(true)}}>Toggle Find</button>
            <button onClick={() => {growState ? setGrowState(false): setGrowState(true)}}>Toggle Grow</button>
            <button onClick={() => {joinState ? setJoinState(false): setJoinState(true)}}>Toggle Join</button>
            <button onClick={() => {linkState ? setLinkState(false): setLinkState(true)}}>Toggle Link</button>
            <button onClick={() => {pruneState ? setPruneState(false): setPruneState(true)}}>Toggle Prune</button>
            <button onClick={() => {rebuildChainState ? setRebuildChainState(false): setRebuildChainState(true)}}>Toggle Rebuild Chain</button>
            <button onClick={() => {sequenceState ? setSequenceState(false): setSequenceState(true)}}>Toggle Sequence</button>
            <button onClick={() => {rebuildBasesState ? setRebuildBasesState(false): setRebuildBasesState(true)}}>Toggle Rebuild Bases</button>

        </div> */}

        <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>
                <MoorhenContainer forwardControls={forwardControls} setMoorhenDimensions={() => {
                    return [1200, 800];
                }} viewOnly={false}/>

            </MoorhenContextProvider>

        </div>

    )

}