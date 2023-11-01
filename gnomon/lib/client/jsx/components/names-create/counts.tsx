import React from "react";
import { makeStyles } from '@material-ui/core/styles';



export interface Count {
    name: string
    description: string
    value: number
    hideAtZero?: boolean
}


const useStyles = makeStyles((theme) => ({
    container: {
        display: "block",
    },
    count: {
        fontStyle: "italic",
        fontWeight: "bold",
    },
    separator: {
        fontWeight: "bold",
    }
}));


const Counts = ({ counts, className, separator, onClick }: {
    counts: Count[],
    className?: string,
    separator?: string,
    onClick?: (countName: string) => void,
}) => {
    const classes = useStyles()

    counts = counts.filter(count => !(count.hideAtZero && count.value == 0))

    return (
        <div
            className={`${classes.container} ${className ? className : ""}`}
        >
            {counts.map((count, idx) =>
                <React.Fragment key={count.name}>
                    <div
                        className={`${classes.count} count count-${count.name}`}
                        onClick={() => onClick && onClick(count.name)}
                    >
                        {`${count.value} ${count.description}`}
                    </div>
                    {separator != undefined
                        && idx < counts.length - 1
                        && counts.length > 1
                        && <div className="separator">{separator}</div>}
                </React.Fragment>
            )}
        </div>
    )
}


export default Counts