import React, { useRef, useEffect } from "react";
import { makeStyles, Theme } from '@material-ui/core/styles';
import InputBase from "@material-ui/core/InputBase";



interface StyleProps {
    minWidth: string
    maxWidth: string
}


const useStyles = makeStyles<Theme, StyleProps>((theme) => ({
    inputContainer: props => ({
        // width: "0",
        minWidth: props.minWidth,
        maxWidth: props.maxWidth,
        "& input::-webkit-outer-spin-button, & input::-webkit-inner-spin-button": {
            display: "none",
            margin: "0",
        },
    }),
    ruler: {
        position: "absolute",
        visibility: "hidden",
        height: "auto",
        width: "auto",
        whiteSpace: "nowrap",
    },
}));


// TODO: make sizing more exact
const AutosizeTextInput = ({ value, onChange, className, type, inputMode = "none", id, placeholder, minWidth = "2em", maxWidth = "unset", inputProps }: {
    value?: string,
    onChange: (event: React.ChangeEvent<HTMLInputElement>) => void,
    className?: string,
    type: string,
    inputMode?: "none" | "text" | "search" | "tel" | "url" | "email" | "numeric" | "decimal"
    id?: string,
    placeholder?: string,
    minWidth?: string,
    maxWidth?: string,
    inputProps: object,
}) => {
    const hasValue = value != undefined
    const rulerRef = useRef<HTMLDivElement>(null)
    const inputContainerRef = useRef<HTMLElement>(null)

    const classes = useStyles({
        minWidth: minWidth,
        maxWidth: maxWidth,
    })

    // handle setting width
    useEffect(() => {
        let inputValue = value
        if (!inputValue) {
            inputValue = placeholder ? placeholder : ""
        }
        
        setInputWidth(inputValue)
    }, [value])

    const setInputWidth = (inputValue: string) => {
        const rulerEl = rulerRef.current
        const inputContainerEl = inputContainerRef.current

        if (rulerEl && inputContainerEl) {
            rulerEl.innerHTML = inputValue
            const widthPx = String(rulerEl.clientWidth)

            inputContainerEl.style.width = `${widthPx}px`
        }
    }
    
    const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        onChange(event)
    }

    return <React.Fragment>
        <InputBase
            value={hasValue ? String(value) : ""}
            onChange={handleInputChange}
            className={`${classes.inputContainer} ${className || ""}`}
            type={type}
            inputMode={inputMode}
            id={id}
            placeholder={placeholder}
            inputProps={inputProps}
            ref={inputContainerRef}
        />
        <div
            className={classes.ruler}
            ref={rulerRef}
        />
    </React.Fragment>
}

export default AutosizeTextInput