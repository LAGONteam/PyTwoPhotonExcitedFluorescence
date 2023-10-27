import talk_to_elliptec_devices as motor
import coherent_laser as laser

FINISHED = False

while FINISHED == False:
    print("list of commands \n",
          "CHANGING POWER LASER: \t",
          "rot_angle (min 0, max 359), example => rot_0 for 0°, rot_47 for 47°... \n",
          "DISPLACING OPO mirror: \t",
          "lin_distance (min 0, max 59), example => lin_0 for home, line_2 for 2mm ... \n",
          "OPEN OPTICAL SHUTTER: \t",
          "shut_state, (open or close only), example => shut_close... \n",
          "OPEN LASER SHUTTER: \t",
          "lasershut_state (open or close only), example => laser_shut_open ... \n",
          "LASER WAVELENGTH TUNING: \t",
          "wave_value (680-1080, only integers), example => wave_951 ... \n",
          "EXIT: \t",
          "exit"
          )
    entry = input("Please send command \n")
    if entry=="exit":
        FINISHED == True
        break
    index = entry.find("_")
    command = entry[:index]
    command_value = entry[index+1:]

    if command == "rot":
        motor.RotationMount().spin_to_position(float(command_value))
    elif command == "lin":
        motor.Linear_stage().displacement(command_value)
    elif command == "shut":
        if command_value == "open":
            motor.Optical_Shutter().open()
        elif command_value == "close":
            motor.Optical_Shutter().close()
    elif command == "lasershut":
        if command_value == 'open':
            laser.Chameleon().openShutter()
        elif command_value == 'close':
            laser.Chameleon().closeShutterBlocking()
    elif command == "wave":
        laser.Chameleon().setWavelengthBlocking(int(command_value))