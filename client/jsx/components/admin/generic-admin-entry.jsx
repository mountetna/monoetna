import * as React from 'react';

export default class GenericAdminEntry extends React.Component{

  constructor(props){

    super(props);

    this['state'] = {

      'editShow': false,
      'editActive': false
    };
  }

  activateEntryEdit(){

    this.setState({ 'editActive': true });
  }

  resetEntry(){

  }

  deactivateEntryEdit(){

    this.setState({ 'editActive': false });
  }

  saveEntry(){

  }

  deleteEntry(){

    
  }

  renderEditControlGroup(){

    // General control group properties.
    var controlBoxProps = {

      'className': 'admin-control-box'
    };

    var editBtnProps = {

      'className': 'admin-entry-edit-btn',
      'title': 'Edit the entry.',
      'onClick': this['activateEntryEdit'].bind(this)
    };

    var resetBtnProps = {

      'className': 'admin-entry-edit-btn',
      'title': 'Reset the entry.',
      'onClick': this['resetEntry'].bind(this)
    };

    var cancelBtnProps = {

      'className': 'admin-entry-edit-btn',
      'title': 'Cancel the entry edit.',
      'onClick': this['deactivateEntryEdit'].bind(this)
    };

    var saveBtnProps = {

      'className': 'admin-entry-edit-btn',
      'title': 'Save the entry edit.',
      'onClick': this['saveEntry'].bind(this)
    };

    var deleteBtnProps = {

      'className': 'admin-entry-delete-btn',
      'title': 'Delete the entry.',
      'onClick': this['deleteEntry'].bind(this)
    };

    // Display settings
    if(this['state']['editShow']){

      controlBoxProps['style'] = { 'visibility': 'visible' };
    }
    else{

      controlBoxProps['style'] = { 'visibility': 'hidden' };
    }

    if(this['state']['editActive']){

      editBtnProps['style'] = { 'display': 'none' };
      resetBtnProps['style'] = { 'display': 'none' };
      cancelBtnProps['style'] = { 'display': 'inline-block' };
      saveBtnProps['style'] = { 'display': 'inline-block' };
    }
    else{

      editBtnProps['style'] = { 'display': 'inline-block' };
      resetBtnProps['style'] = { 'display': 'inline-block' };
      cancelBtnProps['style'] = { 'display': 'none' };
      saveBtnProps['style'] = { 'display': 'none' };
    }

    return (

      <div { ...controlBoxProps }>
                
        <button { ...editBtnProps }>

          <span className='glyphicon glyphicon-pencil'></span>
        </button>
        <button { ...resetBtnProps }>
        
          <span className='glyphicon glyphicon-refresh'></span>
        </button>
        <button { ...cancelBtnProps }>
        
          <span className='glyphicon glyphicon-remove'></span>
        </button>
        <button { ...saveBtnProps }>
        
          <span className='glyphicon glyphicon-ok'></span>
        </button>
         <button { ...deleteBtnProps }>
        
          <span className='glyphicon glyphicon-remove-circle'></span>
        </button>
      </div>
    );
  }

  showControlGroup(event){

    this.setState({ 'editShow': true });
  }

  hideControlGroup(event){

    // Only hide if we are NOT in edit mode.
    if(!this['state']['editActive']) this.setState({ 'editShow': false });
  }

  // This function is overwritten by the inheriting class.
  render(){

    return (<div></div>);
  }
}