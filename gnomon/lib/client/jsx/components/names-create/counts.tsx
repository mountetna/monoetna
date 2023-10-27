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
    },
}));


const Counts = ({ counts, className, separator }: { counts: Count[], className?: string, separator?: string }) => {
    const classes = useStyles()

    return (
        <div className={`${classes.container} ${className ? className : ""}`}>
            {counts.map((count, idx) =>
                <React.Fragment key={count.name}>
                    <div className={`${classes.count} count-${count.name}`}>
                        {
                            count.hideAtZero && count.value == 0
                                ? undefined
                                : `${count.value} ${count.description}`
                        }
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