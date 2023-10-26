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
    },
    count: {
    },
}));


const Counts = ({ counts, className }: { counts: Count[], className?: string }) => {
    const classes = useStyles()

    return (
        <div className={classes.container + className ? ` ${className}` : ""}>
            {counts.map(count =>
                <div className={`${classes.count} count-${count.name}`} key={count.name}>
                    {
                        count.hideAtZero && count.value == 0
                            ? undefined
                            : `${count.value} ${count.description}`
                    }
                </div>
            )}
        </div>
    )
}


export default Counts