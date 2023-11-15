import React, { useState } from "react";
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import InsertDriveFileOutlinedIcon from '@material-ui/icons/InsertDriveFileOutlined';

import ToolbarButtonWithPopper from "./toolbar-button-with-popper";
import { FileFormat, exportDataToBlob } from "../../utils/export";
import { Status } from "../../utils/models";



const ExportButton = ({ small, data, buttonText }: { small: boolean, data: Array<any>, buttonText: string }) => {
    const [open, setOpen] = useState<boolean>(false)
    const [exportStatus, setExportStatus] = useState<Status>("idle")

    const formats: FileFormat[] = ["xlsx", "csv", "tsv"]

    const handleClose = () => {
        setOpen(false);
    };

    const handleExportFile = async (fileFormat: FileFormat) => {
        setExportStatus("inProgress")

        try {
            const blob = await exportDataToBlob(data, fileFormat)
            const filename = `names-${Date.now()}.${fileFormat}`

            const a = document.createElement("a")
            a.setAttribute("href", window.URL.createObjectURL(blob))
            a.setAttribute("download", filename);
            a.click()
            a.remove()

            setExportStatus("idle")
        } catch (err) {
            console.error(`Error export names file: ${err}`)
            setExportStatus("error")
        } finally {
            handleClose()
        }
    }

    return (
        <ToolbarButtonWithPopper
            text={buttonText}
            iconComponent={<InsertDriveFileOutlinedIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            popperComponent={
                <MenuList autoFocusItem={open} id="export-file-formats">
                    {
                        formats.map((format) =>
                            <MenuItem
                                onClick={() => handleExportFile(format)}
                                key={format}
                                disableRipple
                            >
                                {format}
                            </MenuItem>
                        )
                    }
                </MenuList>
            }
            popperId="export-file-formats"
            onClickOrPopperChange={(open: boolean) => setOpen(open)}
            popperOpen={open}
        />
    )
};

export default ExportButton;