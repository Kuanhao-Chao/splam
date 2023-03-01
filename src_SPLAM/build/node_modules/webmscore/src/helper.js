
import LibMscore from '../webmscore.lib.js'
import { IS_NODE, getSelfURL } from './utils.js'

const moduleOptions = IS_NODE
    ? {
        locateFile(path) {
            const { join } = require('path')
            return join(__dirname, path)
        },
        getPreloadedPackage(remotePackageName) {
            const buf = require('fs').readFileSync(remotePackageName).buffer
            return buf
        }
    }
    : {
        locateFile(path) {
            // %INJECTION_HINT_0%
            // fix loading the preload pack in Browsers and WebWorkers
            const prefix = typeof MSCORE_SCRIPT_URL == 'string'
                ? MSCORE_SCRIPT_URL  // to use like an environment variable
                : getSelfURL()
            return new URL(path, prefix).href
        }
    }

/** @type {Record<string, any>} */
let Module = moduleOptions
/** @type {Promise<any>} */
const ModulePromise = LibMscore(moduleOptions)
export { Module }

/**
 * get the pointer to a js string, as utf-8 encoded char*
 * @param {string} str 
 * @returns {number}
 */
export const getStrPtr = (str) => {
    const maxSize = str.length * 4 + 1
    const buf = Module._malloc(maxSize)
    Module.stringToUTF8(str, buf, maxSize)
    return buf
}

/**
 * get the pointer to a TypedArray, as char*
 * @typedef {Int8Array | Int16Array | Int32Array | Uint8Array | Uint16Array | Uint32Array | Uint8ClampedArray | Float32Array | Float64Array} TypedArray
 * @param {TypedArray} data 
 * @returns {number}
 */
export const getTypedArrayPtr = (data) => {
    const size = data.length * data.BYTES_PER_ELEMENT
    const buf = Module._malloc(size)
    Module.HEAPU8.set(data, buf)
    return buf
}

export class WasmRes {
    /**
     * Read responses from the wasm module
     * @param {number} ptr char* pointer to the responses data
     */
    constructor(ptr) {
        /** @type {number} */
        this._ptr = ptr
        /** @type {number} */
        this._size = WasmRes._getUint32(this._sizePtr)
        this._checkRet()
    }

    /**
     * pointer to the error code
     * @private
     */
    get _retCodePtr() {
        return this._ptr
    }

    /**
     * pointer to the data size 
     * @private
     */
    get _sizePtr() {
        return this._retCodePtr + 4
    }

    /**
     * pointer to the data contents 
     * @private
     */
    get _dataPtr() {
        return this._sizePtr + 4
    }

    /**
     * @private 
     * throw error if not ok
     */
    _checkRet() {
        const retCode = WasmRes._getUint32(this._retCodePtr)
        if (retCode !== WasmError.CODE_OK) {
            // read the error message from data
            const retMsg = this.text()
            this.free()
            throw new WasmError(retCode, retMsg)
        }
    }

    /**
     * Read the data contents as Uint8Array
     * @returns {Uint8Array}
     */
    data() {
        return new Uint8Array( // make a copy
            Module.HEAPU8.subarray(this._dataPtr, this._dataPtr + this._size)
        )
    }

    /**
     * Read the data contents as UTF-8 string
     * @returns {string}
     */
    text() {
        return Module.UTF8ToString(this._dataPtr)
    }

    /**
     * Read the data contents as number
     * @returns {number}
     */
    number() {
        return WasmRes._getUint32(this._dataPtr)
    }

    free() {
        return freePtr(this._ptr)
    }

    /**
     * @private 
     * @param {number} ptr 
     * @returns {number}
     */
    static _getUint32(ptr) {
        const sizeData = new DataView(
            new Uint8Array(  // make a copy
                Module.HEAPU8.subarray(ptr, ptr + 4)
            ).buffer
        )
        return sizeData.getUint32(0, true)
    }

    /**
     * @private 
     * @param {number} ptr 
     * @param {keyof WasmRes} method
     */
    static _readAndFree(ptr, method) {
        const res = new WasmRes(ptr)
        const s = res[method]()
        res.free()
        return s
    }

    /**
     * read wasm responses as Uint8Array  
     * @param {number} ptr 
     * @returns {Uint8Array}
     */
    static readData(ptr) {
        return WasmRes._readAndFree(ptr, 'data')
    }

    /**
     * read wasm responses as UTF-8 string 
     * @param {number} ptr 
     * @returns {string}
     */
    static readText(ptr) {
        return WasmRes._readAndFree(ptr, 'text')
    }

    /**
     * read wasm responses as number
     * @param {number} ptr 
     * @returns {number}
     */
    static readNum(ptr) {
        return WasmRes._readAndFree(ptr, 'number')
    }
}

/**
 * free a pointer
 * @param {number} bufPtr 
 */
export const freePtr = (bufPtr) => {
    Module._free(bufPtr)
}

/**
 * this promise is resolved when the runtime is fully initialized
 */
export const RuntimeInitialized = new Promise((resolve) => {
    ModulePromise.then((_Module) => {
        Module = _Module
        Module.ccall('init')  // init libmscore
        resolve(undefined)
    })
})

export class WasmError extends Error {
    /**
     * @param {number} errorCode 
     * @param {string} msg
     */
    constructor(errorCode, msg) {
        super()
        this.name = 'WasmError'
        this.errorCode = errorCode
        this.errorName = msg
        this.message = `WebMscore Err${this.errorName}`
    }
}

/** @type {0} */
WasmError.CODE_OK = 0
