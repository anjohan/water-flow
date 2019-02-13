function tell {
    funnet=$(($(squeue | grep -E "aji[0-9]" | wc -l) + 0))
}

function vent {
    sleep 10
    tell
    while [ ! $funnet -eq 0 ]
    do
        sleep 10
        tell
    done
    echo "Ferdig!"
}
