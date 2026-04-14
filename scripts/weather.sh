#!/usr/bin/env bash
for CITY in "A_Coruna" "Carballo"; do
  echo "--- ${CITY//_/ } ---"
  DATA=$(curl -s "wttr.in/${CITY}?format=j1")

  echo "Morning:   $(echo $DATA | jq -r '.weather[0].hourly[3] | "\(.tempC)°C, \(.lang_es[0].value // .weatherDesc[0].value)"')"
  echo "Midday:    $(echo $DATA | jq -r '.weather[0].hourly[4] | "\(.tempC)°C, \(.lang_es[0].value // .weatherDesc[0].value)"')"
  echo "Afternoon: $(echo $DATA | jq -r '.weather[0].hourly[6] | "\(.tempC)°C, \(.lang_es[0].value // .weatherDesc[0].value)"')"

  if [ "$CITY" == "A_Coruna" ]; then
    HOURLY=$(echo $DATA | jq -r '.weather[0].hourly[4]')
    echo "🌊 Swell:  $(echo $HOURLY | jq -r '.swellHeight_m')m ($(echo $HOURLY | jq -r '.swellPeriod_secs')s) dir $(echo $HOURLY | jq -r '.swellDir16Point')"
    echo "💨 Wind:   $(echo $HOURLY | jq -r '.windspeedKmph')km/h dir $(echo $HOURLY | jq -r '.winddir16Point')"
  fi
  echo ""
done
