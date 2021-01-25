import React, { Component } from 'react';

var ColorPicker = React.createClass({
  render: function() {
    return <div className="color_picker">
             <span className="label">{ this.props.label }</span>
             <input type='text'/>
           </div>
  },

  componentDidMount: function () {
    var picker = "input[type=text]"

    var self = this;
    var node = $(React.findDOMNode(this));
    node.find(picker).spectrum({
      color: this.props.defaultValue,
      showAlpha: true,
      change: function(color) {
        self.props.onChange(color);
      }
    });
  }
});

module.exports = ColorPicker;
