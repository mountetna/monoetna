import sass from 'rollup-plugin-sass'
import babel from "rollup-plugin-babel";
import commonjs from "rollup-plugin-commonjs";
import resolve from "rollup-plugin-node-resolve";
import replace from "rollup-plugin-replace";
import peerDepsExternal from "rollup-plugin-peer-deps-external";

import packageJson from './package.json'

const NODE_ENV = process.env.NODE_ENV || "development";

export default {
  input: "./src/index.js",
  output: {
    file: packageJson.main,
    format: "cjs",
    exports: "named",
    sourcemap: NODE_ENV === "development",
    strict: false
  },
  plugins: [
    peerDepsExternal(),
    commonjs({
        include: 'node_modules/**',
    }),
    babel({
      exclude: "node_modules/**"
    }),
    replace({
      "process.env.NODE_ENV": JSON.stringify(NODE_ENV)
    }),
    resolve(),
    sass({ insert: true })
  ],
  external: [
      'react',
      'react-redux',
      'react-dom',
      'redux'
  ]
};