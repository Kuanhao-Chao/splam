# Changelog

All notable changes to this project will be documented in this file.

### Unreleased

### To be added

* Stream audio file exporting
* Python API
* Pure WebAssembly API

## 1.2.1 - 2023-01-23

### Added

* Export `setLogLevel` on the WebMscore web worker

### Fixed

* Static class fields cause error in babel-transpiled js

## 1.2.0 - 2023-01-23

### Added

* Throw `WasmError` when webmscore encounters score processing error

### Changed

* The log level can now be set before `WebMscore.load(...)`, overrides the default behaviour which turns off logs in `WebMscore.load(...)`.

### Fixed

* For score files created in MuseScore versions older than v3.6.0, instruments no longer all sound like piano in the exported audio.

## 1.1.0 - 2023-01-15

### Added

* Allowing to set log level

```js
const score = await WebMscore.load(...) // set log level to `Off` by default
await WebMscore.setLogLevel(2)          // set log level to `Debug`
```

### Changed

* Logs are turned off by default

* Incorporate [ccache](https://ccache.dev/) in the build script, could significantly speed up recompilation if `ccache` is installed

## 1.0.0 - 2023-01-13

### BREAKING CHANGE

* webmscore is now built based on MuseScore 4.0!

### Changed

* Although there's no public API change comparing to `v0.22.0`, there might be some unintended internal changes or incompatibilities.

## 0.22.0 - 2022-12-31

### Changed

* Build over the `offscreen` [Qt Platform Plugin](https://doc.qt.io/qt-5/qpa.html) so that the wasm binary size has been significantly reduced
    * `webmscore.lib.wasm`:<br>
       10.47 MB -> 7.42 MB
    * `webmscore.lib.mem.wasm`:<br>
       4.97 MB -> 3.92 MB

## 0.21.0 - 2020-03-02

### BREAKING CHANGE

* [`SynthRes.chunk`](https://github.com/LibreScore/webmscore/blob/web/web-public/schemas.ts#L239) becomes **non-interleaved** float32 PCM data

```
SynthRes.chunk = Float32Array[ channelA 512 frames, channelB 512 frames ]
```

## 0.20.4 - 2020-03-01

### Fixed

* Webpack import

## 0.20.2 - 2021-02-08

### Changed

* Improve `synthAudioBatch` performance

## 0.20.0 - 2021-02-08

### Added

* Synthesize audio frames in bulk (n * 512 frames)

```js
const fn = await score.synthAudioBatch(starttime, batchSize)
for (...) {
    const resArr = await fn(cancel?)
}
```

## 0.19.2 - 2021-01-25

### Changed

* Emit detailed debug info (function names) in a separate `.symbols` file

This saves 2.45 MB (23%) uncompressed, 0.42 MB (13%) brotli'd for the generated wasm file, compared to `v0.19.1`.

<details>

<summary>How to look up function names in the .symbols file?</summary>

```log
RuntimeError: function signature mismatch
    at wasm-function[1078]:0x315ac
    at wasm-function[28467]:0x99f067
    ...
```

`webmscore.lib.symbols`

```symbols
1078:Ms::Element::abbox() const
...
28467:Ms::savePng(Ms::Score*, QIODevice*, int, bool, bool)
```

</details>

## 0.19.1 - 2021-01-25

### Added

* Support loading Guitar Pro files into webmscore

```js
await WebMscore.load(format: 'gtp' | 'gp3' | 'gp4' | 'gp5' | 'gpx' | 'gp' | 'ptb', data, fonts, doLayout)
```

### Changed

* Generate detailed debug info

<details>

```log
// Before
RuntimeError: function signature mismatch
    at <anonymous>:wasm-function[1117]:0x359ec
    at <anonymous>:wasm-function[28469]:0x99e040
    ...

// After
RuntimeError: function signature mismatch
    at Ms::Element::abbox() const (wasm-function[1117]:0x359ec)
    at Ms::savePng(Ms::Score*, QIODevice*, int, bool, bool) (wasm-function[28469]:0x99e040)
    ...
```

</details>

## 0.19.0 - 2021-01-23

### Added

* Support loading MusicXML (`musicxml` or compressed `mxl`) and MIDI files into webmscore

```js
await WebMscore.load(format: 'musicxml' | 'mxl' | 'midi', data, fonts, doLayout)
```

## 0.18.0 - 2021-01-20

### Added

* Support MSCZ/MSCX internal file format version 3.02

### Changed

* Merge [upstream v3.6.0](https://github.com/musescore/MuseScore/releases/tag/v3.6) into webmscore, see the changelogs there

## 0.17.1 - 2021-01-12

### Changed

* The `fonts` parameter in `WebMscore.load` can now be passed as type of `Promisable<Uint8Array[]>`.

This allows webmscore to initiate fonts and wasm assets simultaneously.

## 0.17.0 - 2021-01-12

### Changed

* For Webpack build, asset files (`webmscore.lib.*`) will be served via the jsDelivr CDN by default. (importing `webmscore` resolves to `webmscore/webmscore.cdn.mjs`)  
So those files will no longer appear in your `dist` directory, and the workaround is no longer needed.

## 0.16.0 - 2021-01-05

### Added

* Support raw PCM audio data export (for future internal use) 

### Changed

* The original error stack of worker errors is reported in the property `originalStack`
* The workaround for Webpack 4 is no longer enabled by default. To enable it, run `cd node_modules/webmscore && npm run install:webpack` after package installation.

### Removed

* The legacy browser support is removed

## 0.15.2 - 2020-12-31

### Fixed

* [Class fields compatibility issues](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Classes/Public_class_fields#Browser_compatibility)

## 0.15.1 - 2020-12-28

### Fixed

* Remove the previous SoundFont file when `setSoundFont` is called multiple times

## 0.15.0 - 2020-12-28

### Added

* Expose `setSoundFont` as a static method for Node.js  

This allows loading soundfonts in advance, before any score is loaded.

## 0.14.3 - 2020-12-28

### Fixed

* SoundFont related bug ([`2b07ee9`](https://github.com/LibreScore/webmscore/commit/2b07ee9095d243e831afc8128914af73a80228d7))

## 0.14.2 - 2020-12-22

### Changed 

* Data files are in `.wasm` file extension for Webpack. This makes web servers compress those files automatically.

## 0.14.1 - 2020-12-22

### Changed 

* WebMscore web worker no longer extends the native `Worker` class

## 0.14.0 - 2020-12-18

### Added

* Support legacy browsers that have no WebAssembly support

## 0.13.4 - 2020-12-18

### Fixed

* ES5 compilation

## 0.13.3 - 2020-12-18

### Added

* ES5 Compatibility (Webpack only)

## 0.13.2 - 2020-12-16

### Fixed

* Handle exceptions sent from the webmscore web worker

## 0.13.1 - 2020-12-16

### Fixed

* Fixed an issue that webmscore catches all (and irrelevant) uncaught exceptions in Node.js

## 0.13.0 - 2020-12-14

### Added

* Allow to set custom `stdout` and `stderr`  
(Node.js exclusive feature)

See the example in [web-public/src/index.js](https://github.com/LibreScore/webmscore/blob/e7a9f3bed0059842e7ca758f66e75ee8b6ccbd1a/web-public/src/index.js#L34-L67)

## 0.12.1 - 2020-12-14

### Fixed

* Throw file loading error instead of failing silently

## 0.12.0 - 2020-12-14

### Added

* Support MP3 export 

```js
await score.saveAudio('mp3')
```

## 0.11.0 - 2020-12-13

### Added

* The `instrumentName` field of parts in Metadata JSON

## 0.10.3 - 2020-08-10

### Added

* Commonjs Compatibility

## 0.10.0 - 2020-08-01

### BREAKING CHANGE

* The return value of `synthAudio`'s iterator function has changed, see [`interface SynthRes` in schemas.ts](https://github.com/LibreScore/webmscore/blob/web/web-public/schemas.ts#L213)

## 0.9.1 - 2020-07-31

no production code change

## 0.9.0 - 2020-07-30

### BREAKING CHANGE

* The `.destory()` method is now having the `soft` parameter (default: `true`) 
    * `true`: destroy the score instance only, or
    * `false`: destroy the whole WebMscore context

    > To retrieve the default `.destory()` behavior of webmscore in WebWorker prior to v0.9.0, set `soft` to `false`

### Fixed

* Don't transfer the ownership of SoundFont data to the webworker context

## 0.8.3 - 2020-07-30

### Added

* Webpack Compatibility

### Changed 

* Use WOFF2 fonts and LZ4 compression to reduce the data file size (significantly)

## 0.8.1 - 2020-07-28

### Added

* TypeScript declaration (`.d.ts`) files

## 0.8.0 - 2020-07-28

### Added

* The `pageSize` (in pixels) field in position JSON (`measurePositions()` or `segmentPositions()`)

## 0.7.0 - 2020-05-31

### Added

* Save part score as MSCZ/MSCX file

* Boost Mode (set the `doLayout` parameter in `WebMscore.load` to `false`) if you only need score metadata or midi file

### Fixed

* Fixed the runtime error in `processSynth`

## 0.6.0 - 2020-05-29

### Added

* Generate audio files in WAV, OGG, or FLAC formats

* Synthesize raw audio frames, can be used in the Web Audio API 

> A soudfont (sf2/sf3) file must be set using `setSoundFont`

### Changed 

* CJK fonts are no longer bundled inside webmscore. If you would like to export sheet PDF/images with those characters, pass an array of font file (ttf/ttc/otf/otc/woff/etc.) data (Uint8Array) to the `fonts` parameter in `WebMscore.load`

### Fixed

* Always load scores in page mode

## 0.5.0 - 2020-05-13

### Added

* Generate excerpts (if no excerpts found) from parts using `await score.generateExcerpts()`

## 0.4.0 - 2020-05-13

### Changed

* The `name` param is replaced by `filetype` (`"mscz"` or `"mscx"`) in `WebMscore.load`

## 0.3.0 - 2020-05-13

### Added 

* Export individual parts (or `excerpt`s, linked parts). Set the excerpt range of exported files using `await score.setExcerptId(n)`.

* Excerpt information is available in score metadata

## 0.2.1 - 2020-05-12

> Changelog of early versions are omitted

## 0.2.0 - 2020-05-11

## 0.1.0 - 2020-05-10

