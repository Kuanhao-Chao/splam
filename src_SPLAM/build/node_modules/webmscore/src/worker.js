
// The main entry point inside a Web Worker

import WebMscore from './index.js'

/** @type {WebMscore} */
let score

/**
 * @typedef {{ id: number; method: Exclude<keyof import('./index').default, 'scoreptr' | 'excerptId'> | 'load' | 'ready' | 'setLogLevel'; params: any[]; }} RPCReq
 * @typedef {{ id: number; result?: any; error?: any; }} RPCRes
 * @param {number} id 
 * @param {any} result 
 * @param {Transferable[]} transfer
 */
const rpcRes = (id, result, transfer = undefined) => {
    /** @type {RPCRes} */
    const res = {
        id,
        result,
    }
    self.postMessage(res, transfer)
}

/**
 * @param {number} id 
 * @param {Error} err 
 */
const rpcErr = (id, err) => {
    /** @type {RPCRes} */
    const res = {
        id,
        error: {
            name: err.name,
            message: err.message,
            stack: err.stack,
        },
    }
    self.postMessage(res)
}

/**
 * @typedef {import('../schemas').SynthRes | Uint8Array | undefined} Res
 * @param {Res | Res[]} obj 
 * @returns {Transferable[] | undefined}
 */
const getTransferable = (obj) => {
    if (!obj) return
    if (Array.isArray(obj)) {
        return obj.reduce((p, c) => p.concat(getTransferable(c)), []).filter(Boolean)
    } else if (obj instanceof Uint8Array) {
        return [obj.buffer]
    } else if (obj.chunk instanceof Uint8Array) {
        return [obj.chunk.buffer]
    }
}

self.onmessage = async (e) => {
    /** @type {RPCReq} */
    const req = e.data  // JSON-RPC
    const { id, method, params } = req

    try {
        switch (method) {
            case 'ready':
                await WebMscore.ready
                rpcRes(id, 'done')
                break

            case 'load':
                await WebMscore.ready
                score = await WebMscore.load.apply(undefined, params)
                rpcRes(id, 'done')
                break;

            case 'setLogLevel':
                await WebMscore.setLogLevel.apply(undefined, params)
                rpcRes(id, 'done')
                break

            default:
                if (!score) { rpcErr(id, new Error('Score not loaded')) }
                const result = await score[method].apply(score, params)
                rpcRes(id, result, getTransferable(result))
        }
    } catch (err) {
        rpcErr(id, err)
    }
}

