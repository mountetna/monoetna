var config = require('./webpack.config')('development');
var path = require('path');
var fs = require('fs');

config.entry = ['./lib/client/jsx/dev-app.tsx', './lib/client/scss/application.scss']
config.output = null;

config.devServer = {
  contentBase: path.join(__dirname, 'public'),
  compress: true,
  port: 9000,
  host: '0.0.0.0',
};

fs.writeFileSync(path.join(__dirname, "public", "index.html"), fs.readFileSync(path.join(__dirname, "devserver.html")))

module.exports = config;