import * as React from 'react';

export default class ListHead extends React.Component{

  constructor(){

    super();
  }

  selectFile(){

    /*
     * We are using a button to surragate the file input so we may have 
     * a custom browse button.
     */
    document.getElementById('file-selector').click();
  }

  fileSelected(event){

    if(event === undefined) return;
    var fileSelector = event.target;
    var file = fileSelector.files[0];
    this['props'].fileSelected(file);
  }

  render(){

    var fileSelector = {

      'id': 'file-selector',
      'className': 'file-selector',
      'type': 'file',
      'name': 'upload-file',
      'onChange': this['fileSelected'].bind(this)
    };

    return (  

      <thead>

        <tr id='list-head-group'>

          <th id='list-type-column' className='list-head-title'>

            { 'type ' }
            <div className='list-column-head-arrow-group'>

              <span className='glyphicon glyphicon-triangle-bottom'></span>
            </div>
          </th>
          <th id='list-name-column' className='list-head-title'>

            { 'file name ' }
            <div className='list-column-head-arrow-group'>

              <span className='glyphicon glyphicon-triangle-bottom'></span>
            </div>
          </th>
          <th id='list-project-column' className='list-head-title'>

            { 'project ' }
            <div className='list-column-head-arrow-group'>

              <span className='glyphicon glyphicon-triangle-bottom'></span>
            </div>
          </th>
          <th id='list-size-column' className='list-head-title'>

            { 'size ' }
            <div className='list-column-head-arrow-group'>

              <span className='glyphicon glyphicon-triangle-bottom'></span>
            </div>
          </th>
          <th id='list-control-column' className='list-head-title'>
          
            <input { ...fileSelector } />
            <button id='file-select-btn' onClick={ this['selectFile'].bind(this) }>

              <span className='glyphicon glyphicon-plus white-glyphicon'></span>
              { ' ADD FILE' }
            </button>
          </th>
        </tr>
      </thead>
    );
  }
}