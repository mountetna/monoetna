import React from 'react';

// https://reactjs.org/docs/error-boundaries.html
// For #768, because just resetting the selection doesn't
//    happen early enough in the CodeMirror lifecycle...?
export default class ErrorBoundary extends React.Component {
  props: any;
  state: any;
  constructor(props: any) {
    super(props);
    this.state = {hasError: false};
  }

  static getDerivedStateFromError(error: any) {
    // Update state so the next render will show the fallback UI.
    return {hasError: true};
  }
  componentDidCatch(error: any, errorInfo: any) {
    console.error(error);
  }
  render() {
    // We'll still show the Component
    return this.props.children;
  }
}
