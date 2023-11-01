import React from "react";
import { makeStyles, Theme } from '@material-ui/core/styles';
import InputBase from "@material-ui/core/InputBase";


interface StyleProps {
    widthEm: string
    minWidth: string
    maxWidth: string
}


const useStyles = makeStyles<Theme, StyleProps>((theme) => ({
    inputContainer: props => ({
        width: `${props.widthEm}em`,
        minWidth: props.minWidth,
        maxWidth: props.maxWidth,
    }),
}));


// TODO: make it actually sized perfectly or use a 3rd-party lib
const AutosizeInput = ({ value, onChange, className, type, id, placeholder, minWidth = "2em", maxWidth = "8em", inputProps }: {
    value?: string,
    onChange: (event: React.ChangeEvent<HTMLInputElement>) => void,
    className?: string,
    type: string,
    id: string,
    placeholder?: string,
    minWidth?: string,
    maxWidth?: string,
    inputProps: object,
}) => {
    const hasValue = value != undefined

    const classes = useStyles({
        widthEm: String(hasValue ? String(value).length : "0"),
        minWidth: minWidth,
        maxWidth: maxWidth,
    })


    return <InputBase
        value={hasValue ? String(value) : ""}
        onChange={onChange}
        className={`${classes.inputContainer} ${className != undefined ? className : ""}`}
        type={type}
        id={id}
        placeholder={placeholder}
        inputProps={inputProps}
    />
}

export default AutosizeInput