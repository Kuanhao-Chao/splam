
type BoolString = 'true' | 'false'

interface ScorePartData {
    name: string;

    /**
     * MIDI Program Number
     */
    program: number;

    /**
     * @see share/instruments/instruments.xml
     */
    instrumentId: string;

    /**
     * instrument track name
     */
    instrumentName: string;

    lyricCount: number;
    harmonyCount: number;

    hasPitchedStaff: BoolString;
    hasTabStaff: BoolString;
    hasDrumStaff: BoolString;
    isVisible: BoolString;
}

interface ScoreTextFramesData {
    titles: string[];
    subtitles: string[];
    composers: string[];
    poets: string[];
}

interface ScorePageFormat {
    /**
     * page height in mm
     */
    height: number;

    /**
     * page width in mm
     */
    width: number;

    twosided: BoolString;
}

interface ScoreExcerptData {
    /**
     * excerpt id
     */
    id: number;

    /**
     * title of the excerpt
     */
    title: string;

    /**
     * linked parts in the excerpt
     */
    parts: ScorePartData[];
}

/**
 * The score metadata
 */
export interface ScoreMetadata {
    title: string;
    subtitle: string;

    composer: string;

    /** The poet/lyricist of the score */
    poet: string;

    /**
     * The current (lib)mscore version being used
     */
    mscoreVersion: string;

    /**
     * The MSCZ/MSCX file format version
     */
    fileVersion: number;

    /**
     * Number of pages
     */
    pages: number;

    /**
     * Number of measures
     */
    measures: number;

    hasLyrics: BoolString;
    hasHarmonies: BoolString;

    /**
     * @todo explanations
     */
    keysig: number;

    /**
     * ```
     * metaTag("source")
     * ```
     */
    previousSource: string;

    /**
     * @todo explanations
     */
    timesig: string;

    /**
     * The score duration in seconds
     */
    duration: number;

    lyrics: string;

    /**
     * tempo in quarter notes (crochets) per second
     */
    tempo: number;
    /**
     * tempo marker which determines the midi tempo.
     */
    tempoText: string;

    /**
     * excerpts (linked parts) 
     */
    excerpts: ScoreExcerptData[];

    parts: ScorePartData[];

    /**
     * page format of the sheet
     */
    pageFormat: ScorePageFormat;

    /**
     * text frames metadata
     */
    textFramesData: ScoreTextFramesData;
}

export interface PositionElement {
    /**
     * element index
     */
    id: number;

    /**
     * The x coordinate (the top-left corner of the page to the top-left corner of the element)
     */
    x: number;
    /**
     * The y coordinate (the top-left corner of the page to the top-left corner of the element)
     */
    y: number;

    /**
     * The width of the element
     */
    sx: number;
    /**
     * The height of the element
     */
    sy: number;

    /**
     * The page index (zero-based) where the measure or segment presents
     */
    page: number;
}

interface PositionEvent {
    /**
     * The element index corresponding to the envent
     */
    elid: number;

    /**
     * The time position (in ms) of the element (measure or segment) in the exported audio
     */
    position: number;
}

/**
 * The position information of measures or segments
 */
export interface Positions {
    /**
     * The position information of each element (measure or segment) on the sheet,  
     * in pixels of the exported SVG/PNG/PDF file
     * 
     * the `space` property in space.jsonp
     */
    elements: PositionElement[];

    /**
     * The time position/offset (in ms) of each element in the exported audio
     * 
     * the `time` property in space.jsonp
     */
    events: PositionEvent[];

    /**
     * The page size of the exported SVG/PNG/PDF file in pixels
     */
    pageSize: {
        height: number;
        width: number;
    };
}

export interface SynthRes {
    /**
     * Has the value `false` if the iterator is able to produce the next chunk
     */
    done: boolean;

    /**
     * The chunk's starting time (seconds since the score starts)
     */
    startTime: number;

    /**
     * The chunk's end time (seconds since the score starts)  
     * Also the current play time of the synth iterator function
     */
    endTime: number;

    /**
     * The data chunk of audio frames (non-interleaved float32 PCM, 512 frames)  
     * 44100 Hz (44.1 kHz), 0.0116 s (`512 / 44100`)
     */
    chunk: Uint8Array;
}

export type InputFileFormat =
    | 'mscz'             // compressed MuseScore native format
    | 'mscx'             // uncompressed MuseScore native format
    | 'mxl'              // compressed MusicXML
    | 'musicxml' | 'xml' // uncompressed MusicXML
    | 'midi'
    | 'kar'              // unofficial extension of midi
    | 'gtp' | 'gp3' | 'gp4' | 'gp5' | 'gpx' | 'gp' | 'ptb' // Guitar Pro
