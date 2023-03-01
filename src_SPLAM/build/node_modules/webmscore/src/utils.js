
export const IS_NODE = typeof process === 'object' && typeof process.versions === 'object' && typeof process.versions.node === 'string'

export const getSelfURL = () => {
    let url = import.meta.url  // transforms to "" in the generated bundle
    if (!url) {
        if (typeof document !== 'undefined') {
            url = document.currentScript && document.currentScript.src || document.baseURI
        } else if (typeof location !== 'undefined') {
            url = location.href
        }
    }
    return url
}

export const shimDom = () => {
    const getGlobalThis = () => {
        if (typeof self !== 'undefined') { return self }
        if (typeof window !== 'undefined') { return window }
        if (typeof global !== 'undefined') { return global }
        throw new Error('unable to locate global object')
    }

    const globalthis = getGlobalThis()
    globalthis.window = globalthis.window || {
        addEventListener() { },
        location: new URL("file:///"),
        encodeURIComponent,
    }
    globalthis.navigator = globalthis.navigator || {}
}
