
// @ts-check

import {
    Module,
    RuntimeInitialized,
    getStrPtr,
    getTypedArrayPtr,
    WasmRes,
    freePtr,
} from './helper.js'


/** @see WebMscore.hasSoundfont */
let _hasSoundfont = false
/**
 * Don't turn off logs if already set log level before `WebMscore.load(...)` is called
 * @see WebMscore.setLogLevel
 */
let _hasLogLevelSet = false

class WebMscore {

    /**
     * This promise is resolved when the runtime is fully initialized
     * @returns {Promise<void>}
     */
    static get ready() {
        return RuntimeInitialized
    }

    /**
     * The maximum MSCZ/MSCX file format version supported by webmscore 
     * @returns {Promise<number>} e.g. `301`
     */
    static async version() {
        await WebMscore.ready
        return Module.ccall('version', 'number')
    }

    /**
     * Set log level
     * @param {0 | 1 | 2} level - See https://github.com/LibreScore/webmscore/blob/v1.0.0/src/framework/global/thirdparty/haw_logger/logger/log_base.h#L30-L33
     *  - 0: Off
     *  - 1: Normal (`ERRR` or `WARN` or `INFO`)
     *  - 2: Debug  (`DEBG`)
     * @returns {Promise<void>}
     */
    static async setLogLevel(level) {
        _hasLogLevelSet = true
        await WebMscore.ready
        return Module.ccall('setLogLevel', null, ['number'], [level])
    }

    /**
     * Set custom stdout instead of `console.log`  
     * Available before `WebMscore.ready`
     * @private Node.js exclusive
     * @param {(byte: number) => any} write
     */
    static set stdout(write) {
        Module.stdout = write
    }
    /** @private */
    static get stdout() {
        return Module.stdout
    }

    /**
     * Set custom stderr instead of `console.warn`  
     * Available before `WebMscore.ready`
     * @private Node.js exclusive
     * @param {(byte: number) => any} write
     * @example
     * ```js
     * WebMscore['stderr'] = function (byte) {
     *     process.stderr.write(new Uint8Array([byte]))
     * }
     * await WebMscore.ready
     * ```
     */
    static set stderr(write) {
        Module.stderr = write
    }
    /** @private */
    static get stderr() {
        return Module.stderr
    }

    /**
     * Load score data
     * @param {import('../schemas').InputFileFormat} format 
     * @param {Uint8Array} data 
     * @param {Uint8Array[] | Promise<Uint8Array[]>} fonts load extra font files (CJK characters support)
     * @param {boolean} doLayout set to false if you only need the score metadata or the midi file (Super Fast, 3x faster than the musescore software)
     * @returns {Promise<WebMscore>}
     */
    static async load(format, data, fonts = [], doLayout = true) {
        const [_fonts] = await Promise.all([
            fonts,
            WebMscore.ready
        ])

        for (const f of _fonts) {
            await WebMscore.addFont(f)
        }

        const fileformatptr = getStrPtr(format)
        const dataptr = getTypedArrayPtr(data)

        // get the pointer to the MasterScore class instance in C
        const resptr = Module.ccall('load',  // name of C function
            'number',  // return type
            ['number', 'number', 'number', 'boolean'],  // argument types
            [fileformatptr, dataptr, data.byteLength, doLayout]  // arguments
        )
        freePtr(fileformatptr)
        freePtr(dataptr)
        const scoreptr = WasmRes.readNum(resptr)

        if (!_hasLogLevelSet) {
            // turn off logs by default
            await WebMscore.setLogLevel(0);
        }

        const mscore = new WebMscore(scoreptr)
        return mscore
    }

    /**
     * Load (CJK) fonts on demand
     * @private
     * @param {string | Uint8Array} font
     *        * path to the font file in the virtual file system, or
     *        * the font file data
     * @returns {Promise<boolean>} success
     */
    static async addFont(font) {
        if (typeof font !== 'string') {
            const name = '' + Math.random()  // a random name
            // save the font data to the virtual file system
            Module['FS_createDataFile']('/fonts/', name, font, true, true)
            font = '/fonts/' + name
        }

        const fontpathptr = getStrPtr(font)
        const success = Module.ccall('addFont', 'number', ['number'], [fontpathptr])
        freePtr(fontpathptr)
        return !!success
    }

    /**
     * A soundfont file is loaded  
     * @private
     * @type {boolean}
     * @see setSoundFont and saveAudio
     */
    static get hasSoundfont() {
        return _hasSoundfont
    }
    /** @private */
    static set hasSoundfont(value) {
        _hasSoundfont = value
    }

    /**
     * Set the soundfont (sf2/sf3) data  
     * (Audio needs soundfonts)
     * @private
     * @param {Uint8Array} data 
     * @returns {Promise<void>}
     */
    static async setSoundFont(data) {
        if (WebMscore.hasSoundfont) {
            // remove the old soundfont file
            Module['FS_unlink']('/MuseScore_General.sf3')
        }

        // put the soundfont file into the virtual file system
        // side effects: the soundfont is shared across all instances
        Module['FS_createDataFile']('/', 'MuseScore_General.sf3', data, true, true)

        WebMscore.hasSoundfont = true
    }

    /**
     * @hideconstructor use `WebMscore.load`
     * @param {number} scoreptr the pointer to the MasterScore class instance in C++
     */
    constructor(scoreptr) {
        /** @private */
        this.scoreptr = scoreptr

        /** @private */
        this.excerptId = -1
    }

    /**
     * Only save this excerpt (linked parts) of the score  
     * 
     * if no excerpts, generate excerpts from existing instrument parts
     * 
     * @param {number} id  `-1` means the full score 
     */
    async setExcerptId(id) {
        this.excerptId = id
    }

    async getExcerptId() {
        return this.excerptId
    }

    /**
     * Generate excerpts from Parts (only parts that are visible) if no existing excerpts
     * @returns {Promise<void>}
     */
    async generateExcerpts() {
        return Module.ccall('generateExcerpts', null, ['number'], [this.scoreptr])
    }

    /**
     * Get the score title
     * @returns {Promise<string>}
     */
    async title() {
        const dataptr = Module.ccall('title', 'number', ['number'], [this.scoreptr])
        return WasmRes.readText(dataptr)
    }

    /**
     * Get the score title (filename safe, replaced some characters)
     */
    async titleFilenameSafe() {
        const title = await this.title()
        return title.replace(/[\s<>:{}"/\\|?*~.\0\cA-\cZ]+/g, '_')
    }

    /**
     * Get the number of pages in the score (or the excerpt if `excerptId` is set)
     * @returns {Promise<number>}
     */
    async npages() {
        const dataptr = Module.ccall('npages', 'number', ['number', 'number'], [this.scoreptr, this.excerptId])
        return WasmRes.readNum(dataptr)
    }

    /**
     * Get score metadata
     * @returns {Promise<import('../schemas').ScoreMetadata>}
     */
    async metadata() {
        return JSON.parse(await this.saveMetadata())
    }

    /**
     * Get the positions of measures
     * @returns {Promise<import('../schemas').Positions>}
     */
    async measurePositions() {
        return JSON.parse(await this.savePositions(false))
    }

    /**
     * Get the positions of segments
     * @returns {Promise<import('../schemas').Positions>}
     */
    async segmentPositions() {
        return JSON.parse(await this.savePositions(true))
    }

    /**
     * Export score as MusicXML file
     * @returns {Promise<string>} contents of the MusicXML file (plain text)
     */
    async saveXml() {
        const dataptr = Module.ccall('saveXml', 'number', ['number', 'number'], [this.scoreptr, this.excerptId])
        return WasmRes.readText(dataptr)
    }

    /**
     * Export score as compressed MusicXML file
     * @returns {Promise<Uint8Array>}
     */
    async saveMxl() {
        const dataptr = Module.ccall('saveMxl', 'number', ['number', 'number'], [this.scoreptr, this.excerptId])
        return WasmRes.readData(dataptr)
    }

    /**
     * Save part score as MSCZ/MSCX file
     * @param {'mscz' | 'mscx'} format 
     * @returns {Promise<Uint8Array>}
     */
    async saveMsc(format = 'mscz') {
        const dataptr = Module.ccall('saveMsc', 'number', ['number', 'boolean', 'number'], [this.scoreptr, format == 'mscz', this.excerptId])
        return WasmRes.readData(dataptr)
    }

    /**
     * Export score as the SVG file of one page
     * @param {number} pageNumber integer
     * @param {boolean} drawPageBackground 
     * @returns {Promise<string>} contents of the SVG file (plain text)
     */
    async saveSvg(pageNumber = 0, drawPageBackground = false) {
        const dataptr = Module.ccall('saveSvg',
            'number',
            ['number', 'number', 'boolean', 'number'],
            [this.scoreptr, pageNumber, drawPageBackground, this.excerptId]
        )
        return WasmRes.readText(dataptr)
    }

    /**
     * Export score as the PNG file of one page
     * @param {number} pageNumber integer
     * @param {boolean} drawPageBackground 
     * @param {boolean} transparent
     * @returns {Promise<Uint8Array>}
     */
    async savePng(pageNumber = 0, drawPageBackground = false, transparent = true) {
        const dataptr = Module.ccall('savePng',
            'number',
            ['number', 'number', 'boolean', 'boolean', 'number'],
            [this.scoreptr, pageNumber, drawPageBackground, transparent, this.excerptId]
        )
        return WasmRes.readData(dataptr)
    }

    /**
     * Export score as PDF file
     * @returns {Promise<Uint8Array>}
     */
    async savePdf() {
        const dataptr = Module.ccall('savePdf', 'number', ['number', 'number'], [this.scoreptr, this.excerptId])
        return WasmRes.readData(dataptr)
    }

    /**
     * Export score as MIDI file
     * @param {boolean} midiExpandRepeats 
     * @param {boolean} exportRPNs 
     * @returns {Promise<Uint8Array>}
     */
    async saveMidi(midiExpandRepeats = true, exportRPNs = true) {
        const dataptr = Module.ccall('saveMidi',
            'number',
            ['number', 'boolean', 'boolean', 'number'],
            [this.scoreptr, midiExpandRepeats, exportRPNs, this.excerptId]
        )
        return WasmRes.readData(dataptr)
    }

    /**
     * Set the soundfont (sf2/sf3) data
     * @param {Uint8Array} data 
     */
    async setSoundFont(data) {
        return WebMscore.setSoundFont(data)
    }

    /**
     * Export score as audio file (wav/ogg/flac/mp3)
     * @param {'wav' | 'ogg' | 'flac' | 'mp3'} format 
     */
    async saveAudio(format) {
        if (!WebMscore.hasSoundfont) {
            throw new Error('The soundfont is not set.')
        }

        const fileformatptr = getStrPtr(format)
        const dataptr = Module.ccall('saveAudio',
            'number',
            ['number', 'number', 'number'],
            [this.scoreptr, fileformatptr, this.excerptId]
        )
        freePtr(fileformatptr)
        return WasmRes.readData(dataptr)
    }

    /**
     * Synthesize audio frames
     * 
     * `synthAudio` is single instance, i.e. you can't have multiple iterators. If you call `synthAudio` multiple times, it will reset the time offset of all iterators the function returned.
     * 
     * @param {number} starttime The start time offset in seconds
     * @returns {Promise<(cancel?: boolean) => Promise<import('../schemas').SynthRes>>} The iterator function, see `processSynth`
     */
    async synthAudio(starttime) {
        const fn = await this._synthAudio(starttime)
        return (cancel) => {
            return this.processSynth(fn, cancel)
        }
    }

    /**
     * Synthesize audio frames in bulk
     * @param {number} starttime - The start time offset in seconds
     * @param {number} batchSize - max number of result SynthRes' (n * 512 frames)
     * @returns {Promise<(cancel?: boolean) => Promise<import('../schemas').SynthRes[]>>}
     */
    async synthAudioBatch(starttime, batchSize) {
        const fn = await this._synthAudio(starttime)
        return (cancel) => {
            return this.processSynthBatch(fn, batchSize, cancel)
        }
    }

    /**
     * Synthesize audio frames
     * @private
     * @todo GC this iterator function
     * @param {number} starttime The start time offset in seconds
     * @returns {Promise<number>} Pointer to the iterator function
     */
    async _synthAudio(starttime = 0) {
        if (!WebMscore.hasSoundfont) {
            throw new Error('The soundfont is not set.')
        }

        const iteratorFnPtr = Module.ccall('synthAudio',
            'number',
            ['number', 'number', 'number'],
            [this.scoreptr, starttime, this.excerptId]
        )

        const success = iteratorFnPtr !== 0
        if (!success) {
            throw new Error('synthAudio: Internal Error.')
        }

        return iteratorFnPtr
    }

    /**
     * Parse struct SynthRes, then free its memory
     * @private
     * @param {number} resptr - pointer to the SynthRes data
     * @returns {import('../schemas').SynthRes}
     */
    _parseSynthRes(resptr) {
        // struct SynthRes in synthres.h
        const done = Module.getValue(resptr + 0, 'i8')
        const startTime = +Module.getValue(resptr + 4, 'float')  // in seconds
        const endTime = +Module.getValue(resptr + 8, 'float')  // in seconds
        const chunkSize = Module.getValue(resptr + 12, 'i32')
        const chunkPtr = resptr + 16

        const chunk = new Uint8Array(  // make a copy
            Module.HEAPU8.subarray(chunkPtr, chunkPtr + chunkSize)
        )

        freePtr(resptr)

        return {
            done: !!done,
            startTime, // The chunk's start time in seconds
            endTime,   // The current play time in seconds (the chunk's end time)
            chunk,     // The data chunk of audio frames, non-interleaved float32 PCM, 512 frames, 44100 Hz (44.1 kHz), 0.0116 s (512/44100)
        }
    }

    /**
     * @private
     * @param {number} fnptr - pointer to the iterator function
     * @param {boolean} cancel - cancel the audio synthesis worklet 
     * @returns {Promise<import('../schemas').SynthRes>}
     */
    async processSynth(fnptr, cancel = false) {
        const resptr = Module.ccall('processSynth',
            'number',
            ['number', 'boolean'],
            [fnptr, cancel]
        )
        return this._parseSynthRes(resptr)
    }

    /**
     * @private
     * @param {number} fnptr - pointer to the iterator function
     * @param {number} batchSize - see `synthAudioBatch`
     * @param {boolean} cancel - cancel the audio synthesis worklet 
     */
    async processSynthBatch(fnptr, batchSize, cancel = false) {
        const resArrPtr = Module.ccall('processSynthBatch',
            'number',
            ['number', 'number', 'boolean'],
            [fnptr, batchSize, cancel]
        )

        /** @type {import('../schemas').SynthRes[]} */
        const arr = []
        for (let i = 0; i < batchSize; i++) {
            // visit the array of pointers to SynthRes data
            const resptr = Module.getValue(resArrPtr + 4 * i, '*') // 32bit WASM, so one pointer is 4 bytes long
            const r = this._parseSynthRes(resptr)
            arr.push(r)
        }

        freePtr(resArrPtr)
        return arr
    }

    /**
     * Export positions of measures or segments (if `ofSegments` == true) as JSON
     * @param {boolean} ofSegments
     * @also `score.measurePositions()` and `score.segmentPositions()`
     * @returns {Promise<string>}
     */
    async savePositions(ofSegments) {
        const dataptr = Module.ccall('savePositions',
            'number',
            ['number', 'boolean', 'number'],
            [this.scoreptr, ofSegments, this.excerptId]
        )
        return WasmRes.readText(dataptr)
    }

    /**
     * Export score metadata as JSON text
     * @also `score.metadata()`
     * @returns {Promise<string>} contents of the JSON file
     */
    async saveMetadata() {
        const dataptr = Module.ccall('saveMetadata', 'number', ['number'], [this.scoreptr])
        return WasmRes.readText(dataptr)
    }

    /**
     * @param {boolean=} soft (default `true`)
     *                 * `true`  destroy the score instance only, or
     *                 * `false` destroy the whole WebMscore context 
     * @returns {void}
     */
    destroy(soft = true) {
        if (!soft) {
            throw new Error('unimplemented')
        }

        Module.ccall('destroy', 'void', ['number'], [this.scoreptr])
        freePtr(this.scoreptr)
    }

}

export default WebMscore
