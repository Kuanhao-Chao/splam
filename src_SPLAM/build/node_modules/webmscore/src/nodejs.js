
// The main entry point for Node.js

import './nodejs_shim.js'
import WebMscore from './index.js'

WebMscore['default'] = WebMscore  // workaround for commonjs exports
export default WebMscore
